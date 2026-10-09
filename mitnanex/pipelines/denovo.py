"""De novo pipeline: extract mitochondrial reads and assemble the mitogenome.

Port of ``mitnanex_denovo.sh``:

1. length filter + subsample            (seqkit, seqtk)
2. sort reads by length                 (seqkit)
3. all-vs-all overlaps                  (minimap2 ava-ont)
4. cluster + select putative mt reads   (MITNANEX, Python)
5. draft unitigs                        (minimap2 + miniasm)
6. recruit reads from the full dataset  (minimap2 + samtools)
7. quality filter                       (filtlong)
8. assembly                             (Flye --meta)
9. pick the mitogenome contig           (Python)
10. recruit final reads + filter        (minimap2 + samtools + filtlong)
11. final figures (sample reads labelled with the final read set)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .. import plots
from ..assembly import gfa_to_fasta, select_flye_contig
from ..external import run, run_pipeline
from ..extract import ExtractionConfig, MitoReadExtractor, plot_read_level
from ..io import iter_fastq
from .base import Pipeline, default_prefix

log = logging.getLogger("mitnanex")


@dataclass
class DenovoConfig:
    """Options of the de novo pipeline (CLI flags in parentheses).

    Attributes:
        reads: Input FASTQ (``-i``).
        outdir: Output folder (``-o``).
        prefix: Output file prefix (``-r``); default = input file name.
        threads: Threads for external tools (``-t``).
        proportion: Subsampling for the overlap step: fraction (<=1) or number of reads (``-p``).
        min_len / max_len: Read length limits before subsampling (``-m`` / ``-M``).
        minimap_kmer: k-mer size of the all-vs-all overlap search (``-K``).
        extraction: Parameters of the read clustering/selection.
        miniasm_stage: Hidden miniasm ``-S`` option: stop the layout after this stage (``-S``).
        min_mapq: Minimum mapping quality when recruiting reads (``-q``).
        keep_percent: Filtlong ``--keep_percent`` (``-k``).
        flye_mode: Flye read-type flag (``-f``).
        seed: Random seed for subsampling.
    """

    reads: Path
    outdir: Path = Path("mitnanex_results")
    prefix: str | None = None
    threads: int = 4
    proportion: float = 0.4
    min_len: int | None = None
    max_len: int | None = None
    minimap_kmer: int = 15
    extraction: ExtractionConfig = field(default_factory=ExtractionConfig)
    miniasm_stage: int | None = 7
    min_mapq: int = 0
    keep_percent: float = 90
    flye_mode: str = "--nano-hq"
    seed: int = 11


class DenovoPipeline(Pipeline):
    """Reference-free extraction of mitochondrial reads + mitogenome assembly."""

    name = "de novo"
    required_tools = ("seqkit", "seqtk", "minimap2", "miniasm", "samtools", "filtlong", "flye")

    def __init__(self, config: DenovoConfig, timestamped: bool = False, resume: bool = False) -> None:
        super().__init__(config.outdir, config.prefix or default_prefix(config.reads), timestamped, resume)
        self.cfg = config
        self.cfg.reads = Path(config.reads).resolve()
        self.extraction = None

        # All file names in one place.
        self.sample = self.path("_sample.sorted.fastq")
        self.paf = self.path(".paf")
        self.putative_reads = self.path("_putative_mt_reads.fasta")
        self.draft_gfa = self.path("_first_draft_asm.gfa")
        self.draft_fasta = self.path("_first_draft_asm.fasta")
        self.recruit_bam = self.path("_recruit.bam")
        self.recruited = self.path("_collected_reads.fastq")
        self.recruited_filt = self.path("_collected_reads.filtlong.fastq")
        self.flye_dir = self.outdir / f"{self.prefix}_flye_asm"
        self.assembly = self.path("_mt_assembly.fasta")
        self.final_reads = self.path("_mt_reads.fastq")

    def steps(self):
        return [
            ("Filter by length, subsample and sort reads (seqkit, seqtk)", self.subsample, [self.sample]),
            ("All-vs-all read overlaps (minimap2)", self.overlaps, [self.paf]),
            ("Cluster reads and select putative mitochondrial reads (MITNANEX)", self.extract, [self.putative_reads]),
            ("Draft unitigs (miniasm)", self.first_assembly, [self.draft_fasta]),
            ("Recruit reads with the draft unitigs (minimap2)",
             lambda: self.recruit(self.draft_fasta, self.recruited, self.recruited_filt), [self.recruited_filt]),
            ("Assemble with Flye", self.flye, [self.flye_dir / "assembly_info.txt"]),
            ("Select the mitogenome contig", self.select_contig, [self.assembly]),
            ("Recruit final reads with the selected contig",
             lambda: self.recruit(self.assembly, self.path("_final_collected_reads.fastq"), self.final_reads),
             [self.final_reads]),
            ("Final figures", self.final_plots, []),
        ]

    # -- steps -------------------------------------------------------------------
    def subsample(self) -> None:
        """Length filter (seqkit seq) -> random subsample (seqtk) -> sort longest first (seqkit sort)."""
        cfg = self.cfg
        seqkit_seq = ["seqkit", "seq", "--threads", cfg.threads]
        if cfg.min_len:
            seqkit_seq += ["--min-len", cfg.min_len]
        if cfg.max_len:
            seqkit_seq += ["--max-len", cfg.max_len]
        seqkit_seq.append(cfg.reads)
        unsorted = self.path("_sample.fastq")
        run_pipeline([seqkit_seq, ["seqtk", "sample", "-s", cfg.seed, "-", cfg.proportion]],
                     stdout=unsorted, log_file=self.tool_log)
        run(["seqkit", "sort", "--threads", cfg.threads, "--by-length", "--reverse", "-o", self.sample, unsorted],
            log_file=self.tool_log)
        unsorted.unlink()
        self.summary["reads_in_sample"] = sum(1 for _ in iter_fastq(self.sample))
        log.info("%d reads in the sample", self.summary["reads_in_sample"])

    def overlaps(self) -> None:
        """All-vs-all overlaps (``fpa`` filtering is now done in pandas, see overlaps.py)."""
        run(["minimap2", "-x", "ava-ont", "-k", self.cfg.minimap_kmer, "-t", self.cfg.threads, "--dual=yes",
             self.sample, self.sample], stdout=self.paf, log_file=self.tool_log)

    def extract(self) -> None:
        """Run the MITNANEX read clustering/selection."""
        extraction_cfg = self.cfg.extraction
        if self.cfg.min_len and not extraction_cfg.min_overlap_len:
            extraction_cfg.min_overlap_len = self.cfg.min_len  # old: fpa --length-lower $min_len
        result = MitoReadExtractor(extraction_cfg).run(self.sample, self.paf, self.outdir, self.prefix)
        self.extraction = result
        self.summary.update(
            min_cluster_coverage=result.min_coverage,
            n_clusters=len(result.clusters),
            n_selected_clusters=int((result.clusters["status"] == "selected (mt)").sum()),
            putative_mt_reads=len(result.mt_read_ids),
        )

    def first_assembly(self) -> None:
        """Assemble the putative mt reads into unitigs with miniasm, then convert to FASTA."""
        miniasm = ["miniasm"]
        if self.cfg.miniasm_stage:
            miniasm += ["-S", self.cfg.miniasm_stage]
        miniasm += ["-f", self.putative_reads, "-"]
        run_pipeline(
            [["minimap2", "-x", "ava-ont", "-t", self.cfg.threads, "--dual=yes", self.putative_reads, self.putative_reads],
             miniasm],
            stdout=self.draft_gfa, log_file=self.tool_log,
        )
        unitigs = gfa_to_fasta(self.draft_gfa, self.draft_fasta)
        unitigs.to_csv(self.path("_first_draft_asm.tsv"), sep="\t", index=False)
        self.summary["draft_unitigs"] = len(unitigs)
        self.summary["draft_unitigs_total_bp"] = int(unitigs["length"].sum())

    def recruit(self, target: Path, reads_out: Path, filtered_out: Path) -> None:
        """Map *all* input reads to ``target``, keep mapped primary reads, then filtlong."""
        bam = self.recruit_bam
        run_pipeline(
            [["minimap2", "-ax", "map-ont", "--secondary=no", "-t", self.cfg.threads, target, self.cfg.reads],
             # -F 0x904: drop unmapped (0x4), secondary (0x100) and supplementary (0x800) records
             ["samtools", "view", "-b", "-q", self.cfg.min_mapq, "-F", "0x904", "-o", bam, "-"]],
            log_file=self.tool_log,
        )
        run(["samtools", "fastq", bam], stdout=reads_out, log_file=self.tool_log)
        bam.unlink()
        run(["filtlong", "--length_weight", "10", "--keep_percent", self.cfg.keep_percent, reads_out],
            stdout=filtered_out, log_file=self.tool_log)
        n_reads = sum(1 for _ in iter_fastq(filtered_out))
        log.info("Recruited reads: %d after filtlong (%s)", n_reads, filtered_out.name)
        self.summary[f"recruited_reads_{target.stem}"] = n_reads

    def flye(self) -> None:
        """Assemble the recruited reads with Flye (metagenome mode handles uneven coverage)."""
        run(["flye", "-t", self.cfg.threads, "--meta", "--no-alt-contigs", self.cfg.flye_mode,
             self.recruited_filt, "-o", self.flye_dir], log_file=self.tool_log)

    def select_contig(self) -> None:
        best = select_flye_contig(self.flye_dir, self.assembly)
        self.summary.update(
            mt_contig=best["seq_name"], mt_contig_length=int(best["length"]),
            mt_contig_coverage=int(best["coverage"]), mt_contig_circular=best["circular"] == "Y",
            mt_contig_collapsed_from=best["collapsed_from"],
        )

    def final_plots(self) -> None:
        """Compare the clustering output with the final read set in composition space."""
        final_ids = {read_id for read_id, _, _ in iter_fastq(self.final_reads)}
        clustering_ids = set(pd.read_csv(self.path("_read_clusters.tsv"), sep="\t").query("selected")["read_id"])
        labels = pd.Series(
            {read_id: "final mt reads" for read_id in final_ids}
            | {read_id: "clustering only (not recruited)" for read_id in clustering_ids - final_ids}
        )
        palette = {"final mt reads": plots.BLUE, "clustering only (not recruited)": plots.ORANGE,
                   "other reads": plots.GRAY}
        plot_read_level(self.sample, labels, "other reads", self.outdir / "plots", self.prefix, step=6,
                        max_reads=self.cfg.extraction.plot_max_reads, k=self.cfg.extraction.kmer_size,
                        palette=palette, title_suffix="_final",
                        suptitle="Individual reads of the sample, labelled with the final mitochondrial read set")
        self.summary["final_mt_reads"] = len(final_ids)
        log.info("Final assembly: %s", self.assembly)
        log.info("Mitochondrial reads: %s", self.final_reads)
