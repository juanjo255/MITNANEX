"""Helpers to inspect and pick assembly results (miniasm GFA, Flye output)."""

from __future__ import annotations

import io
import logging
import tempfile
from pathlib import Path

import pandas as pd

from .external import capture
from .io import PAF_COLUMNS, read_fasta, read_flye_assembly_info, read_gfa_read_placements, read_gfa_segments, write_fasta
from .kmers import gc_content

log = logging.getLogger(__name__)


def gfa_to_fasta(gfa_file: str | Path, fasta_file: str | Path) -> pd.DataFrame:
    """Convert miniasm unitigs (GFA ``S`` lines) to FASTA and return a summary table.

    Replaces ``gfastats --discover-paths -o``. The returned table has, for
    every unitig, its length, GC content and the number of reads placed on it.
    """
    unitigs = read_gfa_segments(gfa_file)
    if unitigs.empty:
        raise RuntimeError(f"Miniasm produced no unitigs ({gfa_file}). Try a lower --miniasm-stage or more reads.")
    write_fasta(unitigs, fasta_file, line_width=80)

    reads_per_unitig = read_gfa_read_placements(gfa_file).groupby("unitig").size()
    summary = unitigs.assign(
        gc=gc_content(unitigs["seq"]),
        n_reads=unitigs["id"].map(reads_per_unitig).fillna(0).astype(int),
        circular=unitigs["id"].str.endswith("c"),  # miniasm names circular unitigs utgXXXXXXc
    ).drop(columns="seq")
    log.info("Miniasm unitigs:\n%s", summary.to_string(index=False))
    return summary


def tandem_period(seq: str, min_fraction: float = 0.4, min_identity: float = 0.9) -> int | None:
    """Length of the repeated unit if ``seq`` is a tandem concatemer of itself, else ``None``.

    Flye sometimes reports a high-copy circular genome as a circular *dimer*
    (two consecutive copies). The contig is aligned to itself (minimap2, all
    chains); an off-diagonal, same-strand alignment covering >= ``min_fraction``
    of the contig with a constant offset ``d`` means the sequence repeats with
    period ``d``.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".fasta") as handle:
        handle.write(f">contig\n{seq}\n")
        handle.flush()
        text = capture(["minimap2", "-c", "-x", "asm20", "-P", handle.name, handle.name])
    if not text.strip():
        return None
    paf = pd.read_csv(io.StringIO(text), sep="\t", header=None, usecols=range(12), names=PAF_COLUMNS)
    offset = (paf["query_start"] - paf["target_start"]).abs()
    hits = paf[
        (paf["strand"] == "+")
        & (offset > 0.05 * len(seq))                        # not the trivial self-diagonal
        & (paf["aln_len"] >= min_fraction * len(seq))
        & (paf["n_match"] / paf["aln_len"] >= min_identity)
    ]
    if hits.empty:
        return None
    # The smallest offset among long off-diagonal hits is the repeat period.
    return int(offset.loc[hits.index].min())


def collapse_tandem_copies(contigs: pd.DataFrame) -> tuple[pd.DataFrame, int | None]:
    """Keep a single copy of a contig that is a tandem concatemer (see :func:`tandem_period`)."""
    seq = contigs["seq"].iloc[0]
    period = tandem_period(seq)
    if period is None:
        return contigs, None
    log.warning("Contig %s is a tandem concatemer (%d bp = %.1f x %d bp unit); keeping one copy.",
                contigs["id"].iloc[0], len(seq), len(seq) / period, period)
    collapsed = contigs.copy()
    collapsed.loc[collapsed.index[0], "seq"] = seq[:period]
    collapsed["length"] = collapsed["seq"].str.len()
    return collapsed, period


def select_flye_contig(flye_dir: str | Path, output_fasta: str | Path, collapse_concatemers: bool = True) -> pd.Series:
    """Pick the putative mitogenome among Flye contigs and write it to FASTA.

    Rule (as in the original ``select_contig.py``): prefer circular contigs;
    among the candidates take the one with the highest coverage, then the
    longest.

    Fixes a bug of the original: when no contig was circular it took the row
    with *label* 0 (``df["#seq_name"][0]``), i.e. the first line of the file,
    instead of the best contig after sorting.

    If ``collapse_concatemers`` is set, a contig made of tandem copies of the
    mitogenome (a Flye circular dimer) is reduced to one copy.

    Returns:
        The ``assembly_info.txt`` row of the selected contig, with
        ``length`` updated and ``collapsed_from`` set if it was collapsed.
    """
    flye_dir = Path(flye_dir)
    info = read_flye_assembly_info(flye_dir / "assembly_info.txt")
    if info.empty:
        raise RuntimeError("Flye produced no contigs.")
    info = info.sort_values(["coverage", "length"], ascending=False)
    log.info("Flye contigs:\n%s", info.to_string(index=False))

    circular = info[info["circular"] == "Y"]
    best = (circular if not circular.empty else info).iloc[0]
    log.info("Selected contig %s (%s bp, coverage %s, circular=%s)", best["seq_name"], best["length"],
             best["coverage"], best["circular"])

    contigs = read_fasta(flye_dir / "assembly.fasta", ids=[best["seq_name"]])
    best = best.copy()
    best["collapsed_from"] = None
    if collapse_concatemers:
        contigs, period = collapse_tandem_copies(contigs)
        if period:
            best["collapsed_from"], best["length"] = best["length"], period
    write_fasta(contigs, output_fasta, line_width=80)
    return best
