"""Canonical k-mer composition profiles computed with NumPy.

This replaces the Rust extension (``src/utils_rs``) of the original tool.
Encoding a read as an integer array and counting k-mers with ``np.bincount``
is fully vectorised, so it is about as fast as the Rust code without needing a
Rust toolchain / maturin build. Unlike the Rust version it does not crash on
ambiguous bases (``N``): k-mers containing them are skipped.

The profile of a read is the frequency of each *canonical* k-mer (the
lexicographically smaller of a k-mer and its reverse complement), so a read
and its reverse complement get the same profile. For k=3 there are 32
canonical trinucleotides.
"""

from __future__ import annotations

from functools import cached_property
from itertools import product

import numpy as np
import pandas as pd

_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Lookup table byte -> code (A=0, C=1, G=2, T=3, anything else=4).
_ENCODE = np.full(256, 4, dtype=np.uint8)
for _code, _base in enumerate(_BASES):
    _ENCODE[ord(_base)] = _ENCODE[ord(_base.lower())] = _code


def reverse_complement(seq: str) -> str:
    """Reverse complement of an upper-case DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seqs: pd.Series) -> pd.Series:
    """GC fraction of each sequence (vectorised with pandas string methods)."""
    upper = seqs.str.upper()
    gc = upper.str.count("G") + upper.str.count("C")
    return gc / seqs.str.len().clip(lower=1)


class KmerProfiler:
    """Compute normalised canonical k-mer frequency profiles.

    Args:
        k: k-mer size. Small k (3) is robust to Nanopore errors and gives a
            compact feature space (32 dimensions).
    """

    def __init__(self, k: int = 3) -> None:
        if k < 1:
            raise ValueError("k must be >= 1")
        self.k = k
        # Weight of each position inside the k-mer: k-mer index = sum(code * 4**(k-1-i)).
        self._weights = 4 ** np.arange(k - 1, -1, -1)

    @cached_property
    def all_kmers(self) -> list[str]:
        """All 4**k k-mers in lexicographic order (index == integer encoding)."""
        return ["".join(p) for p in product(_BASES, repeat=self.k)]

    @cached_property
    def canonical_kmers(self) -> list[str]:
        """Sorted list of canonical k-mers (the profile columns)."""
        return sorted({min(kmer, reverse_complement(kmer)) for kmer in self.all_kmers})

    @cached_property
    def _to_canonical(self) -> np.ndarray:
        """Array mapping each k-mer integer code to its canonical column index."""
        column = {kmer: i for i, kmer in enumerate(self.canonical_kmers)}
        return np.array([column[min(kmer, reverse_complement(kmer))] for kmer in self.all_kmers])

    def counts(self, seq: str) -> np.ndarray:
        """Raw canonical k-mer counts of one sequence."""
        codes = _ENCODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
        n_kmers = len(codes) - self.k + 1
        if n_kmers <= 0:
            return np.zeros(len(self.canonical_kmers))

        # Sliding windows of length k (a view, no copy) -> one integer per k-mer.
        windows = np.lib.stride_tricks.sliding_window_view(codes, self.k)
        valid = (windows < 4).all(axis=1)          # skip k-mers with N/ambiguous bases
        kmer_codes = windows[valid] @ self._weights

        forward_counts = np.bincount(kmer_codes, minlength=4**self.k)
        # Fold each k-mer onto its canonical representative.
        return np.bincount(self._to_canonical, weights=forward_counts, minlength=len(self.canonical_kmers))

    def profile(self, seq: str) -> np.ndarray:
        """Canonical k-mer frequencies of one sequence (sums to 1)."""
        counts = self.counts(seq)
        total = counts.sum()
        return counts / total if total else counts

    def profile_many(self, seqs: pd.Series) -> pd.DataFrame:
        """Profiles of many sequences; rows keep the index of ``seqs``.

        Returns:
            DataFrame (n_sequences x n_canonical_kmers) of frequencies.
        """
        matrix = np.vstack([self.profile(s) for s in seqs]) if len(seqs) else np.empty((0, len(self.canonical_kmers)))
        return pd.DataFrame(matrix, index=seqs.index, columns=self.canonical_kmers)

    def pooled_profiles(self, seqs: pd.Series, groups: pd.Series) -> pd.DataFrame:
        """One profile per group from the summed k-mer counts of its sequences.

        Summing counts before normalising weights every read by its length, so
        a group's profile is the composition of all its bases together.

        Args:
            seqs: Sequences.
            groups: Group label of each sequence (same index as ``seqs``).
        """
        counts = pd.DataFrame(
            np.vstack([self.counts(s) for s in seqs]) if len(seqs) else np.empty((0, len(self.canonical_kmers))),
            index=seqs.index, columns=self.canonical_kmers,
        )
        summed = counts.groupby(groups).sum()
        return summed.div(summed.sum(axis=1).replace(0, 1), axis=0)
