"""Readers and writers for the file formats used by MITNANEX.

All readers return :class:`pandas.DataFrame` objects so the rest of the tool
can rely on vectorised pandas operations (filtering, merging, grouping)
instead of hand-written loops.

Formats covered:
    * FASTQ / FASTA (optionally gzip-compressed, detected from the extension)
    * PAF  (minimap2 pairwise alignment format)
    * GFA1 (miniasm assembly graph; only segment lines are needed)
    * Flye ``assembly_info.txt``
"""

from __future__ import annotations

import gzip
from itertools import islice
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import numpy as np
import pandas as pd

# Column names of the 12 mandatory PAF fields (see minimap2 man page).
PAF_COLUMNS = [
    "query_id",       # 0  query sequence name
    "query_len",      # 1  query sequence length
    "query_start",    # 2  query start (0-based)
    "query_end",      # 3  query end (0-based, exclusive)
    "strand",         # 4  '+' same strand, '-' opposite strand
    "target_id",      # 5  target sequence name
    "target_len",     # 6  target sequence length
    "target_start",   # 7  target start on original strand
    "target_end",     # 8  target end on original strand
    "n_match",        # 9  number of matching bases
    "aln_len",        # 10 alignment block length (including gaps)
    "mapq",           # 11 mapping quality (0-255; 255 = missing)
]

# Compact dtypes keep memory low: all-vs-all PAF files can have millions of rows.
_PAF_DTYPES = {
    "query_id": "string",
    "query_len": np.int32,
    "query_start": np.int32,
    "query_end": np.int32,
    "strand": "category",
    "target_id": "string",
    "target_len": np.int32,
    "target_start": np.int32,
    "target_end": np.int32,
    "n_match": np.int32,
    "aln_len": np.int32,
    "mapq": np.int16,
}


def open_text(path: str | Path) -> TextIO:
    """Open a plain or gzip-compressed text file for reading."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# NOTE on performance: sequence files have few but very long lines.
# ``pd.read_csv`` (used as a line splitter) was benchmarked to be ~12x slower
# than streaming lines with Python's buffered reader (3.1 s vs 0.25 s for
# 160 Mb of reads), so FASTA/FASTQ/GFA are parsed with plain Python and only
# the parsed records are put into a DataFrame. Tabular formats (PAF, VCF,
# Flye's assembly_info) are read with pandas.


def iter_fastq(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(id, sequence, quality)`` for each record of a 4-line FASTQ."""
    with open_text(path) as handle:
        while True:
            record = list(islice(handle, 4))  # one FASTQ record = 4 lines
            if not record:
                return
            if len(record) < 4 or not record[0].startswith("@"):
                raise ValueError(f"{path}: truncated or malformed FASTQ record near {record[0][:50]!r}")
            read_id = record[0][1:].split(None, 1)[0]  # ID = header up to first whitespace
            yield read_id, record[1].rstrip("\n"), record[3].rstrip("\n")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(id, sequence)`` for each record of a (multi-line) FASTA."""
    read_id, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if read_id is not None:
                    yield read_id, "".join(chunks)
                read_id, chunks = line[1:].split(None, 1)[0], []
            else:
                chunks.append(line.strip())
    if read_id is not None:
        yield read_id, "".join(chunks)


def _records_to_frame(records: Iterable[tuple], columns: list[str], ids: Iterable[str] | None) -> pd.DataFrame:
    """Build a DataFrame from parsed records, optionally keeping only ``ids``."""
    if ids is not None:
        wanted = set(ids)
        records = (r for r in records if r[0] in wanted)
    frame = pd.DataFrame.from_records(list(records), columns=columns)
    frame["length"] = frame["seq"].str.len()
    return frame


def read_fastq(path: str | Path, ids: Iterable[str] | None = None, with_quality: bool = False) -> pd.DataFrame:
    """Load a FASTQ file into a DataFrame.

    Args:
        path: FASTQ file (``.gz`` allowed).
        ids: If given, only records whose ID is in this collection are kept
            (the others are never stored, which keeps memory low).
        with_quality: Also return the quality strings (needed to write FASTQ).

    Returns:
        DataFrame with columns ``id``, ``seq``, (``qual``), ``length``.
    """
    frame = _records_to_frame(iter_fastq(path), ["id", "seq", "qual"], ids)
    return frame if with_quality else frame.drop(columns="qual")


def read_fasta(path: str | Path, ids: Iterable[str] | None = None) -> pd.DataFrame:
    """Load a (multi-line) FASTA file into a DataFrame with ``id``, ``seq``, ``length``."""
    return _records_to_frame(iter_fasta(path), ["id", "seq"], ids)


def read_sequences(path: str | Path, ids: Iterable[str] | None = None) -> pd.DataFrame:
    """Read FASTA or FASTQ, deciding from the first character of the file."""
    with open_text(path) as handle:
        first_char = handle.read(1)
    if first_char == "@":
        return read_fastq(path, ids=ids)
    if first_char == ">":
        return read_fasta(path, ids=ids)
    raise ValueError(f"{path}: not a FASTA/FASTQ file")


def write_fasta(records: pd.DataFrame, path: str | Path, line_width: int | None = None) -> Path:
    """Write a DataFrame with ``id`` and ``seq`` columns as FASTA.

    Args:
        records: Sequences to write.
        path: Output file.
        line_width: Wrap sequences at this width (``None`` = one line per sequence).
    """
    seqs = records["seq"]
    if line_width:
        # Insert a newline every ``line_width`` characters.
        seqs = seqs.str.findall(f".{{1,{line_width}}}").str.join("\n")
    text = (">" + records["id"].astype(str) + "\n" + seqs).str.cat(sep="\n")
    Path(path).write_text(text + "\n" if text else "")
    return Path(path)


def write_fastq(records: pd.DataFrame, path: str | Path) -> Path:
    """Write a DataFrame with ``id``, ``seq`` and ``qual`` columns as FASTQ."""
    text = ("@" + records["id"].astype(str) + "\n" + records["seq"] + "\n+\n" + records["qual"]).str.cat(sep="\n")
    Path(path).write_text(text + "\n" if text else "")
    return Path(path)


def read_paf(path: str | Path) -> pd.DataFrame:
    """Read the 12 mandatory columns of a PAF file.

    Optional SAM-like tags (``tp:A:P``, ``cm:i:..``) vary in number between
    lines; ``usecols`` makes pandas ignore them.
    """
    paf = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=range(len(PAF_COLUMNS)),
        names=PAF_COLUMNS,
        dtype=_PAF_DTYPES,
        engine="c",
    )
    return paf


def read_gfa_segments(path: str | Path) -> pd.DataFrame:
    """Return the segment (``S``) lines of a GFA1 file as ``id``/``seq``/``length``.

    Miniasm writes unitig sequences in the ``S`` lines (when run with ``-f``),
    so converting the graph to FASTA is just selecting these rows.
    """
    with open_text(path) as handle:
        segments = [line.split("\t", 3)[1:3] for line in handle if line.startswith("S\t")]
    gfa = pd.DataFrame(segments, columns=["id", "seq"])
    gfa["length"] = gfa["seq"].str.len()
    return gfa


def read_gfa_read_placements(path: str | Path) -> pd.DataFrame:
    """Return miniasm ``a`` lines (``unitig``, ``read_id``): the reads that build each unitig."""
    # Format: a <unitig> <offset> <read:start-end> <strand> <length>
    with open_text(path) as handle:
        rows = [line.split("\t")[1:4:2] for line in handle if line.startswith("a\t")]
    placements = pd.DataFrame(rows, columns=["unitig", "read"])
    placements["read_id"] = placements["read"].str.rsplit(":", n=1).str[0]
    return placements.drop(columns="read")


def read_flye_assembly_info(path: str | Path) -> pd.DataFrame:
    """Read Flye's ``assembly_info.txt`` with clean column names."""
    info = pd.read_csv(path, sep="\t")
    return info.rename(
        columns={
            "#seq_name": "seq_name",
            "cov.": "coverage",
            "circ.": "circular",
            "repeat": "repeat",
            "mult.": "multiplicity",
            "alt_group": "alt_group",
            "graph_path": "graph_path",
        }
    )
