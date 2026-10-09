"""De novo extraction of putative mitochondrial reads (the core of MITNANEX).

Input: a set of reads (FASTQ) and their all-vs-all overlaps (PAF).
Steps:

1. filter overlaps (identity, internal matches) -> :mod:`mitnanex.overlaps`
2. greedy read clustering                       -> :mod:`mitnanex.clustering`
3. keep high-coverage clusters                  -> :func:`coverage_threshold_by_gap`
4. split them by trinucleotide composition      -> :class:`CompositionSelector`
5. write the reads of the selected clusters, tables and diagnostic plots.

This module replaces the old ``main.py``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import plots
from .clustering import ClusteringResult, ReadClusterer
from .io import read_fastq, read_paf, write_fasta
from .kmers import KmerProfiler, gc_content
from .overlaps import OverlapFilter
from .selection import CompositionSelection, CompositionSelector, coverage_threshold_by_gap

log = logging.getLogger(__name__)

MT_LABEL = "retrieved mt"
OTHER_LABEL = "other reads"


@dataclass
class ExtractionConfig:
    """Parameters of the de novo extraction (defaults follow the original tool).

    Attributes:
        min_identity: Minimum overlap identity to cluster two reads together.
        min_overlap_len: Minimum overlap (alignment block) length; ``None`` = off.
        drop_internal: Discard internal matches (shared repeats).
        min_coverage: Fixed minimum cluster coverage; ``None`` = automatic gap.
        min_clusters: Automatic cut-off must keep more than this many clusters.
        kmer_size: k for the composition profiles.
        profile: Composition of a cluster from ``"members"`` (all its reads
            pooled) or ``"longest"`` (its longest read only, original MITNANEX).
        group_rule: Choice of the mitochondrial K-means group: ``"max"`` (group of
            the highest-coverage cluster) or ``"median"`` (original rule).
        haplotype_mode: One K-means group per cluster (old ``-y`` option): only
            the single best cluster is kept, useful to separate haplotypes.
        plot_max_reads: Reads used for read-level composition plots.
        make_plots: Produce the diagnostic figures.
        seed: Random seed (PCA / K-means / read sampling for plots).
    """

    min_identity: float = 0.6
    min_overlap_len: int | None = None
    drop_internal: bool = True
    min_coverage: int | None = None
    min_clusters: int = 3
    kmer_size: int = 3
    profile: str = "members"
    group_rule: str = "max"
    haplotype_mode: bool = False
    plot_max_reads: int = 6000
    make_plots: bool = True
    seed: int = 0


@dataclass
class ExtractionResult:
    """Everything produced by :class:`MitoReadExtractor`."""

    clustering: ClusteringResult
    min_coverage: int
    composition: CompositionSelection
    mt_read_ids: pd.Index
    clusters: pd.DataFrame
    outputs: dict[str, Path] = field(default_factory=dict)


class MitoReadExtractor:
    """Run the de novo extraction and write its outputs.

    Example:
        >>> extractor = MitoReadExtractor(ExtractionConfig(min_identity=0.6))
        >>> result = extractor.run("sample.sorted.fastq", "sample.paf", "results/", prefix="sample")
    """

    def __init__(self, config: ExtractionConfig | None = None) -> None:
        self.config = config or ExtractionConfig()

    # -- algorithm -----------------------------------------------------------
    def cluster_reads(self, paf: pd.DataFrame) -> ClusteringResult:
        """Filter overlaps and cluster the reads."""
        cfg = self.config
        overlap_filter = OverlapFilter(
            min_identity=cfg.min_identity, min_aln_len=cfg.min_overlap_len, drop_internal=cfg.drop_internal
        )
        kept = overlap_filter.apply(paf)
        log.info("Overlaps: %d in PAF, %d kept after filtering", len(paf), len(kept))
        clustering = ReadClusterer().fit(kept)
        log.info("Reads clustered: %d in %d clusters", len(clustering.membership), len(clustering.clusters))
        return clustering

    def select_clusters(self, clustering: ClusteringResult, reads_file: str | Path) -> tuple[int, CompositionSelection]:
        """Coverage filter followed by the composition-based selection."""
        cfg = self.config
        clusters = clustering.clusters
        if clusters.empty:
            raise RuntimeError("No read clusters were formed: check --min-identity and the PAF file.")

        min_coverage = cfg.min_coverage if cfg.min_coverage else coverage_threshold_by_gap(
            clusters["coverage"], cfg.min_clusters
        )
        candidates = clusters[clusters["coverage"] >= min_coverage]
        log.info("Minimum coverage: %d reads -> %d candidate clusters", min_coverage, len(candidates))

        # Reads describing each candidate cluster: all members or just the longest read.
        if cfg.profile == "members":
            reads_of = clustering.membership[clustering.membership["cluster_id"].isin(candidates["cluster_id"])]
        elif cfg.profile == "longest":
            reads_of = candidates.rename(columns={"repr_read_id": "read_id"})[["read_id", "cluster_id"]]
        else:
            raise ValueError("profile must be 'members' or 'longest'")
        seqs = read_fastq(reads_file, ids=reads_of["read_id"]).set_index("id")["seq"]
        cluster_seqs = reads_of.assign(seq=reads_of["read_id"].map(seqs)).dropna(subset=["seq"])
        n_groups = len(candidates) if cfg.haplotype_mode else 2
        selector = CompositionSelector(k=cfg.kmer_size, n_groups=n_groups, random_state=cfg.seed,
                                       group_rule=cfg.group_rule)
        return min_coverage, selector.fit(candidates, cluster_seqs)

    @staticmethod
    def annotate_clusters(clusters: pd.DataFrame, min_coverage: int, composition: CompositionSelection) -> pd.DataFrame:
        """Add a human-readable ``status`` (and composition columns) to every cluster."""
        table = clusters.merge(
            composition.table.drop(columns=["coverage", "repr_read_id", "repr_read_len"]),
            on="cluster_id", how="left",
        )
        table["status"] = np.select(
            [table["selected"].eq(True), table["coverage"] >= min_coverage],
            ["selected (mt)", "rejected by composition"],
            default="below coverage cut-off",
        )
        return table

    # -- entry point -----------------------------------------------------------
    def run(self, reads_file: str | Path, paf_file: str | Path, outdir: str | Path, prefix: str = "mitnanex") -> ExtractionResult:
        """Run the whole extraction and write outputs into ``outdir``.

        Outputs:
            ``{prefix}_putative_mt_reads.fasta``, ``{prefix}_clusters.tsv``,
            ``{prefix}_read_clusters.tsv`` and figures in ``plots/``.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        clustering = self.cluster_reads(read_paf(paf_file))
        min_coverage, composition = self.select_clusters(clustering, reads_file)
        clusters = self.annotate_clusters(clustering.clusters, min_coverage, composition)

        selected_ids = clusters.loc[clusters["status"] == "selected (mt)", "cluster_id"]
        mt_read_ids = clustering.reads_of(selected_ids)
        log.info("MITNANEX retrieved %d putative mitochondrial reads from %d clusters", len(mt_read_ids), len(selected_ids))

        outputs = {
            "mt_reads": write_fasta(read_fastq(reads_file, ids=mt_read_ids), outdir / f"{prefix}_putative_mt_reads.fasta"),
            "clusters": outdir / f"{prefix}_clusters.tsv",
            "read_clusters": outdir / f"{prefix}_read_clusters.tsv",
        }
        clusters.to_csv(outputs["clusters"], sep="\t", index=False, float_format="%.6g")
        clustering.membership.assign(selected=clustering.membership["read_id"].isin(mt_read_ids)).to_csv(
            outputs["read_clusters"], sep="\t", index=False
        )

        result = ExtractionResult(clustering, min_coverage, composition, mt_read_ids, clusters, outputs)
        if self.config.make_plots:
            outputs.update(self.plot(result, reads_file, outdir / "plots", prefix))
        return result

    # -- figures ---------------------------------------------------------------
    def plot(self, result: ExtractionResult, reads_file: str | Path, plot_dir: Path, prefix: str) -> dict[str, Path]:
        """Draw the diagnostic figures for an extraction."""
        figures = {
            "plot_coverage": plots.plot_cluster_coverage(
                result.clusters, result.min_coverage, plot_dir / f"{prefix}_1_cluster_coverage.png"
            ),
            "plot_candidates": plots.plot_candidate_composition(
                result.composition.table, result.composition.explained_variance,
                plot_dir / f"{prefix}_2_candidate_composition.png",
            ),
        }
        labels = pd.Series(MT_LABEL, index=result.mt_read_ids)
        figures.update(
            plot_read_level(reads_file, labels, OTHER_LABEL, plot_dir, prefix, step=3,
                            max_reads=self.config.plot_max_reads, k=self.config.kmer_size, seed=self.config.seed)
        )
        return figures


def sample_reads_for_plot(reads: pd.DataFrame, labels: pd.Series, background: str, max_reads: int, seed: int) -> pd.DataFrame:
    """Label reads and down-sample them to ``max_reads`` keeping every labelled group visible.

    Labelled reads (e.g. retrieved mt) are kept up to half of the budget; the
    rest of the budget is filled with background reads.
    """
    reads = reads.assign(group=reads["id"].map(labels).fillna(background))
    labelled = reads[reads["group"] != background]
    labelled = labelled.sample(n=min(len(labelled), max_reads // 2), random_state=seed)
    others = reads[reads["group"] == background]
    others = others.sample(n=min(len(others), max_reads - len(labelled)), random_state=seed)
    return pd.concat([labelled, others], ignore_index=True)


def plot_read_level(reads_file: str | Path, labels: pd.Series, background: str, plot_dir: Path, prefix: str,
                    step: int, max_reads: int = 6000, k: int = 3, seed: int = 0,
                    palette: dict[str, str] | None = None, title_suffix: str = "",
                    suptitle: str = "Individual reads of the sample") -> dict[str, Path]:
    """Read-level composition, trinucleotide signature and read-length figures.

    Args:
        reads_file: Reads to describe (FASTQ).
        labels: Read ID -> group name for the reads of interest.
        background: Group name for every other read.
        step: Number used to order the figure file names.
        palette: Group -> color (default: blue for the first label group, gray background).
        title_suffix: Added to output keys/file names (e.g. ``"_final"``).
        suptitle: Title of the read-composition figure.
    """
    reads = sample_reads_for_plot(read_fastq(reads_file), labels, background, max_reads, seed)
    if palette is None:
        palette = {name: color for name, color in zip(labels.unique(), [plots.BLUE, plots.ORANGE, plots.AQUA])}
        palette[background] = plots.GRAY

    profiles = KmerProfiler(k).profile_many(reads["seq"])
    coords, explained = plots.read_composition_pca(profiles, random_state=seed)
    reads = pd.concat([reads, coords], axis=1).assign(gc=gc_content(reads["seq"]))
    reads.attrs["explained_variance"] = explained

    figures = {
        f"plot_read_composition{title_suffix}": plots.plot_read_composition(
            reads, palette, plot_dir / f"{prefix}_{step}_read_composition{title_suffix}.png",
            suptitle=suptitle,
        ),
        f"plot_read_lengths{title_suffix}": plots.plot_read_lengths(
            reads, palette, plot_dir / f"{prefix}_{step + 1}_read_lengths{title_suffix}.png"
        ),
    }
    target = next(iter(palette))
    if (reads["group"] == target).any() and (reads["group"] == background).any():
        figures[f"plot_kmer_signature{title_suffix}"] = plots.plot_kmer_signature(
            profiles, reads["group"], target, background, plot_dir / f"{prefix}_{step + 2}_kmer_signature{title_suffix}.png"
        )
    return figures
