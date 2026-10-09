"""Greedy clustering of reads from their all-vs-all overlaps.

Idea: mitochondrial DNA is present in many copies per cell, so mitochondrial
reads overlap with many other reads and end up in large ("high coverage")
clusters, while most single-copy nuclear reads form small clusters.

The algorithm reproduces the original MITNANEX logic (``assign_cluster.py``):
overlaps are visited in file order and for each read pair

1. neither read is clustered  -> open a new cluster with both reads;
2. only one read is clustered -> add the other read to that cluster;
3. both reads are clustered   -> do nothing (clusters are never merged).

The pairwise loop is inherently sequential, so it is a plain Python loop over
NumPy arrays with a dictionary lookup; all per-cluster statistics are then
computed with a single pandas ``groupby``.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class ClusteringResult:
    """Output of :class:`ReadClusterer`.

    Attributes:
        membership: One row per clustered read: ``read_id``, ``cluster_id``,
            ``read_len``. Rows are in the order reads joined their cluster.
        clusters: One row per cluster: ``cluster_id``, ``coverage`` (number of
            reads), ``repr_read_id`` / ``repr_read_len`` (the longest read,
            used to represent the cluster), sorted by decreasing coverage.
    """

    membership: pd.DataFrame
    clusters: pd.DataFrame

    def reads_of(self, cluster_ids) -> pd.Index:
        """IDs of all reads that belong to any of ``cluster_ids``."""
        mask = self.membership["cluster_id"].isin(cluster_ids)
        return pd.Index(self.membership.loc[mask, "read_id"])


class ReadClusterer:
    """Assign reads to clusters from a filtered PAF table (see module docstring)."""

    def fit(self, paf: pd.DataFrame) -> ClusteringResult:
        """Cluster the reads present in ``paf``.

        Args:
            paf: Overlaps already filtered by :class:`~mitnanex.overlaps.OverlapFilter`.
        """
        read_to_cluster: dict[str, int] = {}  # read ID -> cluster ID
        read_length: dict[str, int] = {}      # read ID -> length (from the PAF)
        n_clusters = 0

        # NumPy arrays make the loop ~10x faster than DataFrame.itertuples().
        columns = [paf[c].to_numpy() for c in ("query_id", "target_id", "query_len", "target_len")]
        for query, target, query_len, target_len in zip(*columns):
            query_cluster = read_to_cluster.get(query)
            target_cluster = read_to_cluster.get(target)

            if query_cluster is None and target_cluster is None:
                # Case 1: open a new cluster containing both reads.
                read_to_cluster[query] = read_to_cluster[target] = n_clusters
                read_length[query], read_length[target] = int(query_len), int(target_len)
                n_clusters += 1
            elif target_cluster is None:
                # Case 2a: query already clustered -> recruit the target read.
                read_to_cluster[target] = query_cluster
                read_length[target] = int(target_len)
            elif query_cluster is None:
                # Case 2b: target already clustered -> recruit the query read.
                read_to_cluster[query] = target_cluster
                read_length[query] = int(query_len)
            # Case 3: both reads already clustered -> nothing to do.

        membership = pd.DataFrame(
            {
                "read_id": list(read_to_cluster.keys()),
                "cluster_id": list(read_to_cluster.values()),
                "read_len": list(read_length[r] for r in read_to_cluster),
            }
        )
        return ClusteringResult(membership=membership, clusters=self.summarise(membership))

    @staticmethod
    def summarise(membership: pd.DataFrame) -> pd.DataFrame:
        """Per-cluster coverage and representative (longest) read."""
        if membership.empty:
            return pd.DataFrame(columns=["cluster_id", "coverage", "repr_read_id", "repr_read_len"])
        grouped = membership.groupby("cluster_id", sort=False)
        # idxmax returns the first longest read, matching the original tie-breaking.
        longest = membership.loc[grouped["read_len"].idxmax()].set_index("cluster_id")
        clusters = pd.DataFrame(
            {
                "coverage": grouped.size(),
                "repr_read_id": longest["read_id"],
                "repr_read_len": longest["read_len"],
            }
        ).reset_index()
        return clusters.sort_values(["coverage", "repr_read_len"], ascending=False, ignore_index=True)
