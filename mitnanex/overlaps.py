"""Filtering of all-vs-all read overlaps (PAF) with vectorised pandas code.

This replaces the external ``fpa drop --internalmatch --length-lower`` call
used by the original shell pipeline (``fpa`` has no osx-arm64 build) and the
per-line identity filter that was done inside the Python loop.

Overlap classification follows miniasm (Li 2016, Bioinformatics 32:2103,
Algorithm 5). An *internal match* is an alignment that stops well before the
ends of both reads -- typically a shared repeat rather than a true overlap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class OverlapFilter:
    """Keep only overlaps that are likely to come from the same locus.

    Attributes:
        min_identity: Minimum ``n_match / aln_len`` (the old ``-s`` option).
        min_aln_len: Minimum alignment block length; ``None`` disables it
            (equivalent to fpa ``--length-lower``).
        drop_internal: Drop internal matches (fpa ``--internalmatch``).
        max_hang: Maximum overhang allowed before calling an internal match.
        int_frac: Overhang must also be below ``int_frac * mapping length``.
        drop_self: Remove self-alignments (a read mapped to itself).
    """

    min_identity: float = 0.6
    min_aln_len: int | None = None
    drop_internal: bool = True
    max_hang: int = 1000
    int_frac: float = 0.8
    drop_self: bool = True

    @staticmethod
    def identity(paf: pd.DataFrame) -> pd.Series:
        """Fraction of matching bases in the alignment block."""
        return paf["n_match"] / paf["aln_len"]

    @staticmethod
    def classify(paf: pd.DataFrame, max_hang: int = 1000, int_frac: float = 0.8) -> pd.Series:
        """Label each overlap as internal / contained / dovetail (miniasm algorithm 5).

        Coordinates of the target are flipped for reverse-strand hits so that
        both reads are compared in the same orientation.
        """
        reverse = (paf["strand"] == "-").to_numpy()
        q_start, q_end, q_len = (paf[c].to_numpy(np.int64) for c in ("query_start", "query_end", "query_len"))
        t_len = paf["target_len"].to_numpy(np.int64)
        # Flip target coordinates when the alignment is on the opposite strand.
        t_start = np.where(reverse, t_len - paf["target_end"].to_numpy(np.int64), paf["target_start"].to_numpy(np.int64))
        t_end = np.where(reverse, t_len - paf["target_start"].to_numpy(np.int64), paf["target_end"].to_numpy(np.int64))

        # Unaligned bases on the left and on the right side of the overlap.
        overhang = np.minimum(q_start, t_start) + np.minimum(q_len - q_end, t_len - t_end)
        map_len = np.maximum(q_end - q_start, t_end - t_start)

        internal = overhang > np.minimum(max_hang, map_len * int_frac)
        q_contained = (q_start <= t_start) & (q_len - q_end <= t_len - t_end)
        t_contained = (q_start >= t_start) & (q_len - q_end >= t_len - t_end)

        # np.select evaluates conditions in order: first match wins.
        labels = np.select(
            [internal, q_contained, t_contained],
            ["internal", "query_contained", "target_contained"],
            default="dovetail",
        )
        return pd.Series(labels, index=paf.index, name="overlap_type")

    def apply(self, paf: pd.DataFrame) -> pd.DataFrame:
        """Return the filtered PAF table (original row order is preserved).

        Row order matters because read clustering is greedy and processes
        overlaps in the order minimap2 wrote them.
        """
        keep = pd.Series(True, index=paf.index)
        if self.drop_self:
            keep &= paf["query_id"] != paf["target_id"]
        if self.min_aln_len:
            keep &= paf["aln_len"] >= self.min_aln_len
        keep &= self.identity(paf) >= self.min_identity
        if self.drop_internal:
            keep &= self.classify(paf, self.max_hang, self.int_frac) != "internal"
        return paf.loc[keep]
