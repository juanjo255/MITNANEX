"""Evaluate a `mitnanex denovo` run on reads simulated with ``tests/simulate.py``.

Read names start with their true origin (``mt_``, ``nuc_``, ``rdna_``, ...),
so precision / recall of each read set can be computed, and the composition
plots can be redrawn with the *true* origin of each read.

Usage:
    python tests/benchmark.py --run-dir results/ --prefix sim --reads sim.fastq
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mitnanex import plots
from mitnanex.extract import plot_read_level
from mitnanex.io import read_fastq, read_sequences


def origin(ids: pd.Series) -> pd.Series:
    """True origin encoded in the simulated read name."""
    return ids.str.split("_").str[0]


def score(read_ids: pd.Series, truth_ids: pd.Series) -> dict[str, object]:
    """Precision / recall of a set of reads predicted as mitochondrial."""
    predicted_mt = origin(read_ids) == "mt"
    true_positives = int(predicted_mt.sum())
    return {
        "reads": len(read_ids),
        "precision": round(true_positives / max(len(read_ids), 1), 3),
        "recall": round(true_positives / max(int((origin(truth_ids) == "mt").sum()), 1), 3),
        "contaminants": origin(read_ids[~predicted_mt]).value_counts().to_dict(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--reads", type=Path, required=True, help="the full simulated FASTQ")
    args = parser.parse_args()

    run = args.run_dir
    sample = run / f"{args.prefix}_sample.sorted.fastq"
    sample_ids = read_fastq(sample)["id"]
    all_ids = read_fastq(args.reads)["id"]

    table = pd.DataFrame({
        "clustering (vs. sample)": score(read_sequences(run / f"{args.prefix}_putative_mt_reads.fasta")["id"], sample_ids),
        "final reads (vs. all reads)": score(read_fastq(run / f"{args.prefix}_mt_reads.fastq")["id"], all_ids),
    }).T
    print(table.to_string())
    table.to_csv(run / f"{args.prefix}_benchmark.tsv", sep="\t")

    # Same read-level figures, coloured by the true origin of each read.
    labels = pd.Series(origin(sample_ids).to_numpy(), index=sample_ids.to_numpy())
    labels = labels[labels != "nuc"].replace({"mt": "mitochondrial", "rdna": "rDNA array",
                                              "te": "transposon", "numt": "NUMT"})
    palette = {"mitochondrial": plots.BLUE, "rDNA array": plots.ORANGE, "transposon": plots.AQUA,
               "nuclear (single copy)": plots.GRAY}
    labels = labels[labels.isin(palette)]  # the few NUMT reads are drawn as nuclear
    plot_read_level(sample, labels, "nuclear (single copy)", run / "plots", args.prefix, step=9,
                    palette=palette, title_suffix="_truth",
                    suptitle="Individual reads of the sample, coloured by their TRUE origin (simulation)")


if __name__ == "__main__":
    main()
