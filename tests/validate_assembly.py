"""Validate a `mitnanex denovo` result against a published mitogenome.

* Assembly: aligned with minimap2 (asm presets); reports how much of the
  reference is covered, alignment identity and contig length.
* Reads: all input reads are mapped to the reference mitogenome. Reads mapped
  with MAPQ >= 20 over >= 50% of their length are taken as "reference mt reads";
  precision / recall of MITNANEX's final read set are computed against them
  (an approximation: NUMT reads also map and strain differences reduce MAPQ).

Usage:
    python tests/validate_assembly.py --run-dir results/ --prefix sample \\
        --reference ref_mito.fasta --reads sample.fastq [--threads 8]
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mitnanex.external import capture
from mitnanex.io import PAF_COLUMNS, iter_fastq, read_fasta


def paf_from_text(text: str) -> pd.DataFrame:
    """Parse minimap2 PAF output (mandatory columns only) with pandas."""
    if not text.strip():
        return pd.DataFrame(columns=PAF_COLUMNS)
    return pd.read_csv(io.StringIO(text), sep="\t", header=None, usecols=range(12), names=PAF_COLUMNS)


def _cs_operations(cs_tags: pd.Series) -> pd.DataFrame:
    """Split minimap2 ``cs`` tags into operations: match (``:``), substitution (``*``), ins (``+``), del (``-``)."""
    ops = cs_tags.str.extractall(r"(?P<op>[:*+\-])(?P<value>[0-9]+|[a-z]+)")
    ops["length"] = np.where(ops["op"] == ":", pd.to_numeric(ops["value"], errors="coerce"), ops["value"].str.len())
    ops.loc[ops["op"] == "*", "length"] = 1  # '*ag' is one substituted base
    return ops


def assembly_stats(assembly: Path, reference: Path) -> dict:
    """Reference coverage, error counts and gap-compressed identity of the assembly.

    Gap-compressed identity counts an indel as one difference whatever its
    length, so a large strain-specific insertion/deletion (e.g. an intron)
    does not hide an otherwise accurate assembly.
    """
    text = capture(["minimap2", "-c", "--cs", "-x", "asm20", "--secondary=no", reference, assembly])
    aln = paf_from_text(text)
    ref_len = int(read_fasta(reference)["length"].iloc[0])
    contig_len = int(read_fasta(assembly)["length"].sum())
    if aln.empty:
        return {"contig_length": contig_len, "reference_length": ref_len, "reference_covered": 0.0}
    # Union of covered intervals on the reference (alignments may overlap at a circular junction).
    covered = set()
    for start, end in zip(aln["target_start"], aln["target_end"]):
        covered.update(range(start, end))

    cs = pd.Series([line.split("cs:Z:")[1].split("\t")[0] for line in text.splitlines()])
    ops = _cs_operations(cs)
    counts = ops.groupby("op")["length"].agg(["size", "sum"])
    get = lambda op, col: int(counts.loc[op, col]) if op in counts.index else 0  # noqa: E731
    matches, subs, ins, dels = get(":", "sum"), get("*", "size"), get("+", "size"), get("-", "size")
    indels = ops[ops["op"].isin(["+", "-"])]
    return {
        "contig_length": contig_len,
        "reference_length": ref_len,
        "reference_covered": round(len(covered) / ref_len, 4),
        "substitutions": subs,
        "insertions": ins,
        "deletions": dels,
        "largest_indel_bp": int(indels["length"].max()) if len(indels) else 0,
        "gap_compressed_identity": round(matches / (matches + subs + ins + dels), 5),
        "n_alignments": len(aln),
    }


def read_stats(run_dir: Path, prefix: str, all_reads: Path, reference: Path, threads: int) -> dict:
    """Precision / recall of the recruited reads vs. reads mapping to the reference.

    Computed before filtlong (``collected``) and after it (``final``): filtlong
    keeps a percentage of *bases* favouring long reads, so many short reads are
    dropped on purpose; recall is therefore given in reads and in bases.
    """
    aln = paf_from_text(capture(["minimap2", "-x", "map-ont", "--secondary=no", "-t", threads, reference, all_reads]))
    aln = aln[(aln["mapq"] >= 20) & ((aln["query_end"] - aln["query_start"]) >= 0.5 * aln["query_len"])]
    reference_mt = aln.drop_duplicates("query_id").set_index("query_id")["query_len"]

    stats = {"reference_mt_reads": len(reference_mt), "reference_mt_bases": int(reference_mt.sum())}
    for label, path in [("collected", run_dir / f"{prefix}_final_collected_reads.fastq"),
                        ("final", run_dir / f"{prefix}_mt_reads.fastq")]:
        if not path.exists():
            continue
        ids = pd.Index([read_id for read_id, _, _ in iter_fastq(path)])
        shared = reference_mt.index.intersection(ids)
        stats[label] = {
            "reads": len(ids),
            "precision": round(len(shared) / max(len(ids), 1), 4),
            "recall_reads": round(len(shared) / max(len(reference_mt), 1), 4),
            "recall_bases": round(reference_mt.loc[shared].sum() / max(reference_mt.sum(), 1), 4),
        }
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--reference", type=Path, required=True, help="published mitogenome (FASTA)")
    parser.add_argument("--reads", type=Path, required=True, help="all input reads")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    report = {
        "assembly": assembly_stats(args.run_dir / f"{args.prefix}_mt_assembly.fasta", args.reference),
        "reads": read_stats(args.run_dir, args.prefix, args.reads, args.reference, args.threads),
    }
    print(json.dumps(report, indent=2))
    (args.run_dir / f"{args.prefix}_validation.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
