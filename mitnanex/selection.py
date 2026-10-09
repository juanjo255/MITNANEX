"""Selection of the mitochondrial read clusters.

Two filters are applied to the clusters produced by
:class:`~mitnanex.clustering.ReadClusterer`:

1. **Coverage** - keep clusters above the largest "gap" in the sorted
   coverage values (mitochondrial clusters are expected to stand out), or
   above a user-provided minimum coverage.
2. **Composition** - each surviving cluster is described by the canonical
   trinucleotide frequencies of *all its reads pooled* (default) or of its
   longest read only (original MITNANEX; the longest read of a mitochondrial
   cluster is often a chimeric/NUMT read longer than the mitogenome, whose
   composition misrepresents the cluster), projected to 2-D with
   PCA and split in two groups with K-means (as in metagenomic binning). The
   group with the highest median normalised coverage is taken as
   mitochondrial.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from .kmers import KmerProfiler

log = logging.getLogger(__name__)


def coverage_threshold_by_gap(coverage: pd.Series, min_clusters: int = 3) -> int:
    """Minimum coverage that separates high-coverage clusters from the rest.

    The coverage values are sorted in decreasing order and the largest drop
    between consecutive values is used as the cut-off, provided that more
    than ``min_clusters`` clusters remain above it (otherwise the next largest
    gap is tried). If no gap qualifies, every cluster is kept.

    Args:
        coverage: Coverage (number of reads) per cluster.
        min_clusters: Number of clusters that must be exceeded above the cut-off.

    Returns:
        The coverage value at the upper side of the selected gap.
    """
    cov = coverage.sort_values(ascending=False, ignore_index=True)
    gaps = cov.diff(periods=-1)                       # cov[i] - cov[i+1]
    n_kept = cov.rank(method="max", ascending=False)  # clusters with coverage >= cov[i]

    candidates = pd.DataFrame({"coverage": cov, "gap": gaps, "n_kept": n_kept}).dropna()
    candidates = candidates[candidates["n_kept"] > min_clusters]
    if candidates.empty:
        log.warning("Too few clusters (%d) to find a coverage gap; keeping all of them.", len(cov))
        return int(cov.min()) if len(cov) else 0
    # Stable sort: among equal gaps the one at higher coverage wins.
    best = candidates.sort_values("gap", ascending=False, kind="stable").iloc[0]
    return int(best["coverage"])


@dataclass
class CompositionSelection:
    """Result of :class:`CompositionSelector`.

    Attributes:
        table: One row per candidate cluster with PCA coordinates
            (``PC1``/``PC2``), K-means group, GC content and a boolean
            ``selected`` column.
        profiles: Canonical k-mer frequencies of each candidate cluster.
        pca: The fitted PCA (re-used to project other reads in the plots).
        selected_group: K-means group chosen as mitochondrial.
    """

    table: pd.DataFrame
    profiles: pd.DataFrame
    pca: PCA | None
    selected_group: int | None = None
    explained_variance: np.ndarray = field(default_factory=lambda: np.zeros(2))


class CompositionSelector:
    """PCA + weighted K-means on k-mer profiles of cluster representatives.

    Args:
        k: k-mer size for the composition profile.
        n_groups: Number of K-means groups (2 = mitochondrial vs. other).
        max_iter: K-means iterations.
        random_state: Seed for reproducibility.
        group_rule: How the mitochondrial K-means group is chosen: ``"max"``
            (group of the cluster with the highest normalised coverage) or
            ``"median"`` (highest median normalised coverage, original MITNANEX).
    """

    GROUP_RULES = ("max", "median")

    def __init__(self, k: int = 3, n_groups: int = 2, max_iter: int = 100, random_state: int = 0,
                 group_rule: str = "max") -> None:
        if group_rule not in self.GROUP_RULES:
            raise ValueError(f"group_rule must be one of {self.GROUP_RULES}")
        self.group_rule = group_rule
        self.profiler = KmerProfiler(k)
        self.n_groups = n_groups
        self.max_iter = max_iter
        self.random_state = random_state

    def fit(self, candidates: pd.DataFrame, cluster_seqs: pd.DataFrame) -> CompositionSelection:
        """Select mitochondrial clusters among ``candidates``.

        Args:
            candidates: Cluster table with ``cluster_id``, ``coverage``,
                ``repr_read_id`` and ``repr_read_len``.
            cluster_seqs: Reads describing each cluster (``cluster_id``, ``seq``):
                all member reads, or only the representative read of each cluster.
        """
        table = candidates.copy()
        # Number of reads per base of the representative read: a length-independent
        # proxy for depth, used both as K-means weight and to pick the mt group.
        table["coverage_norm"] = table["coverage"] / table["repr_read_len"]

        # Pooled composition (and GC) of the reads of each cluster.
        profiles = self.profiler.pooled_profiles(cluster_seqs["seq"], cluster_seqs["cluster_id"])
        profiles = profiles.reindex(table["cluster_id"].to_numpy())
        gc_bases = cluster_seqs["seq"].str.count("[GCgc]").groupby(cluster_seqs["cluster_id"]).sum()
        all_bases = cluster_seqs["seq"].str.len().groupby(cluster_seqs["cluster_id"]).sum()
        table["gc"] = table["cluster_id"].map(gc_bases / all_bases)

        # With fewer than 3 clusters there is nothing to separate: keep them all.
        if len(table) < 3:
            log.warning("Only %d candidate cluster(s); skipping composition clustering.", len(table))
            table[["PC1", "PC2"]] = np.nan
            table["kmeans_group"] = 0
            table["selected"] = True
            return CompositionSelection(table=table, profiles=profiles, pca=None, selected_group=0)

        pca = PCA(n_components=2, random_state=self.random_state)
        table[["PC1", "PC2"]] = pca.fit_transform(profiles.to_numpy())

        kmeans = KMeans(
            n_clusters=min(self.n_groups, len(table)),
            max_iter=self.max_iter,
            init="k-means++",
            n_init=1,
            random_state=self.random_state,
        )
        table["kmeans_group"] = kmeans.fit_predict(table[["PC1", "PC2"]], sample_weight=table["coverage_norm"])

        # Mitochondrial DNA is expected to have the highest depth: pick the group
        # holding the top cluster, or (original rule) the highest group median.
        if self.group_rule == "max":
            selected_group = int(table.loc[table["coverage_norm"].idxmax(), "kmeans_group"])
        else:
            selected_group = int(table.groupby("kmeans_group")["coverage_norm"].median().idxmax())
        table["selected"] = table["kmeans_group"] == selected_group
        log.info(
            "PCA explained variance: %s; selected K-means group %d (%d clusters)",
            np.round(pca.explained_variance_ratio_, 3),
            selected_group,
            table["selected"].sum(),
        )
        return CompositionSelection(
            table=table,
            profiles=profiles,
            pca=pca,
            selected_group=selected_group,
            explained_variance=pca.explained_variance_ratio_,
        )
