"""Figures that make every MITNANEX decision visible.

Every function takes tidy pandas tables, draws one figure and saves it as PNG.
Group labels are free text, so the same functions are used for the tool's
own output ("retrieved mt" vs "other reads") and for benchmarks with known
ground truth.

Color rules: the group of interest (mitochondrial) is always blue and is drawn
on top; the background (nuclear / other) is a recessive gray. When a third
group is needed (e.g. clusters passing coverage but rejected by composition)
it gets orange. Text never takes a series color.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping

import matplotlib

matplotlib.use("Agg")  # file output only; no display needed (works on clusters/servers)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

log = logging.getLogger(__name__)

# --- palette (validated categorical slots + neutral background) ---------------
BLUE = "#2a78d6"      # mitochondrial / selected
ORANGE = "#eb6834"    # second group (e.g. rejected by composition)
AQUA = "#1baf7a"      # third group
RED = "#e34948"       # negative pole of diverging bars
GRAY = "#b4b2ab"      # background: nuclear / filtered out
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e6e5e0"

STYLE = {
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": GRID,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "axes.axisbelow": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "legend.frameon": False,
    "text.color": INK,
}


def _save(fig: plt.Figure, path: str | Path) -> Path:
    """Save and close a figure, creating the folder if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    log.info("Plot written: %s", path)
    return path


def _group_order(groups: pd.Series, palette: Mapping[str, str]) -> list[str]:
    """Groups to draw, background first so the group of interest stays on top."""
    present = set(groups.unique())
    return [g for g in reversed(list(palette)) if g in present]


def _scatter_groups(ax, data: pd.DataFrame, x: str, y: str, group: str, palette: Mapping[str, str],
                    size: pd.Series | float = 14, alpha: float = 0.75) -> None:
    """Scatter plot with one color per group and counts in the legend."""
    for name in _group_order(data[group], palette):
        subset = data[data[group] == name]
        s = size.loc[subset.index] if isinstance(size, pd.Series) else size
        ax.scatter(subset[x], subset[y], s=s, c=palette[name], alpha=alpha,
                   edgecolors="white", linewidths=0.4, label=f"{name} (n={len(subset):,})")
    # Legend in palette order (group of interest first).
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: list(palette).index(labels[i].rsplit(" (", 1)[0]))
    # Legend below the axes so it never hides data points.
    legend = ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper left",
                       bbox_to_anchor=(0, -0.16), ncol=1, borderaxespad=0)
    for handle in legend.legend_handles:  # same marker size in the legend whatever the point sizes
        handle.set_sizes([30])
        handle.set_alpha(1)


# -----------------------------------------------------------------------------
# 1. Coverage filter
# -----------------------------------------------------------------------------
def plot_cluster_coverage(clusters: pd.DataFrame, min_coverage: float, path: str | Path) -> Path:
    """Clusters ranked by coverage, showing the coverage cut-off and the final choice.

    Args:
        clusters: Cluster table with ``coverage`` and ``status`` columns
            (``selected`` / ``rejected by composition`` / ``below coverage``).
        min_coverage: Coverage threshold used by the filter.
    """
    palette = {"selected (mt)": BLUE, "rejected by composition": ORANGE, "below coverage cut-off": GRAY}
    data = clusters.sort_values("coverage", ascending=False, ignore_index=True)
    data["rank"] = np.arange(1, len(data) + 1)

    with plt.rc_context(STYLE):
        fig, (ax_rank, ax_hist) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1.6, 1]})
        _scatter_groups(ax_rank, data, "rank", "coverage", "status", palette, size=16)
        ax_rank.axhline(min_coverage, color=INK_SECONDARY, linewidth=0.9)
        ax_rank.annotate(f"cut-off = {min_coverage:g} reads", (data["rank"].max(), min_coverage),
                         xytext=(0, 4), textcoords="offset points", ha="right", color=INK_SECONDARY)
        ax_rank.set_xscale("log")
        ax_rank.set_yscale("log")
        ax_rank.set(xlabel="Cluster rank (by coverage)", ylabel="Coverage (reads per cluster)",
                    title="Read clusters ranked by coverage")

        bins = np.unique(np.logspace(0, np.log10(max(data["coverage"].max(), 2)), 30).round())
        ax_hist.hist(data["coverage"], bins=bins, color=GRAY, edgecolor="white", linewidth=0.5)
        ax_hist.axvline(min_coverage, color=INK_SECONDARY, linewidth=0.9)
        ax_hist.set_xscale("log")
        ax_hist.set_yscale("log")
        ax_hist.set(xlabel="Coverage (reads per cluster)", ylabel="Number of clusters",
                    title="Coverage distribution")
        return _save(fig, path)


# -----------------------------------------------------------------------------
# 2. Composition of candidate clusters (what the algorithm actually sees)
# -----------------------------------------------------------------------------
def plot_candidate_composition(table: pd.DataFrame, explained_variance, path: str | Path) -> Path:
    """PCA of trinucleotide profiles of candidate clusters + GC vs. coverage.

    Args:
        table: ``CompositionSelection.table`` (PC1, PC2, gc, coverage,
            coverage_norm, selected).
        explained_variance: PCA explained-variance ratio (2 values).
    """
    palette = {"selected (mt)": BLUE, "rejected by composition": ORANGE}
    data = table.assign(group=np.where(table["selected"], "selected (mt)", "rejected by composition"))
    # Marker area proportional to cluster coverage (readable range 20-220 pt^2).
    cov = data["coverage"].astype(float)
    size = 20 + 200 * (cov - cov.min()) / max(cov.max() - cov.min(), 1)

    with plt.rc_context(STYLE):
        fig, (ax_pca, ax_gc) = plt.subplots(1, 2, figsize=(10, 4.2))
        if data["PC1"].notna().all():
            _scatter_groups(ax_pca, data, "PC1", "PC2", "group", palette, size=size)
            ev = np.asarray(explained_variance) * 100
            ax_pca.set(xlabel=f"PC1 ({ev[0]:.1f}% variance)", ylabel=f"PC2 ({ev[1]:.1f}% variance)")
        ax_pca.set_title("Trinucleotide composition (PCA)")
        fig.suptitle("Candidate clusters (composition of their reads; marker size = coverage)",
                     x=0.02, ha="left", fontweight="bold", fontsize=11)

        _scatter_groups(ax_gc, data, "gc", "coverage_norm", "group", palette, size=size)
        ax_gc.set_yscale("log")
        ax_gc.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
        ax_gc.set(xlabel="GC content of the cluster", ylabel="Normalised coverage (reads / bp)",
                  title="GC content vs. normalised coverage")
        return _save(fig, path)


# -----------------------------------------------------------------------------
# 3. Composition of individual reads: retrieved mt vs. nuclear background
# -----------------------------------------------------------------------------
def read_composition_pca(profiles: pd.DataFrame, random_state: int = 0) -> tuple[pd.DataFrame, np.ndarray]:
    """2-D PCA of read k-mer profiles; returns coordinates and explained variance."""
    pca = PCA(n_components=2, random_state=random_state)
    coords = pd.DataFrame(pca.fit_transform(profiles.to_numpy()), index=profiles.index, columns=["PC1", "PC2"])
    return coords, pca.explained_variance_ratio_


def plot_read_composition(reads: pd.DataFrame, palette: Mapping[str, str], path: str | Path,
                          suptitle: str = "Individual reads") -> Path:
    """Scatter of read composition (PCA) and GC content vs. read length.

    Args:
        reads: One row per read with ``PC1``, ``PC2``, ``gc``, ``length`` and
            ``group`` columns (see :func:`read_composition_pca`).
        palette: Ordered mapping group -> color; first entry = group of interest.
        suptitle: Figure title above both panels.
    """
    ev = reads.attrs.get("explained_variance", [np.nan, np.nan])
    with plt.rc_context(STYLE):
        fig, (ax_pca, ax_gc) = plt.subplots(1, 2, figsize=(10.5, 4.4))
        _scatter_groups(ax_pca, reads, "PC1", "PC2", "group", palette, size=9, alpha=0.6)
        ax_pca.set(xlabel=f"PC1 ({ev[0] * 100:.1f}% variance)", ylabel=f"PC2 ({ev[1] * 100:.1f}% variance)",
                   title="Trinucleotide composition (PCA)")
        fig.suptitle(suptitle, x=0.02, ha="left", fontweight="bold", fontsize=11)

        _scatter_groups(ax_gc, reads, "length", "gc", "group", palette, size=9, alpha=0.6)
        ax_gc.set_xscale("log")
        ax_gc.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
        ax_gc.set(xlabel="Read length (bp)", ylabel="GC content", title="GC content vs. read length")
        return _save(fig, path)


def plot_kmer_signature(profiles: pd.DataFrame, groups: pd.Series, target: str, background: str,
                        path: str | Path) -> Path:
    """Which trinucleotides distinguish the target reads from the background.

    Bars show log2(mean frequency in ``target`` / mean frequency in
    ``background``) for every canonical trinucleotide; blue = enriched in the
    target group, red = depleted.
    """
    mean_freq = profiles.groupby(groups).mean()
    log2_ratio = np.log2(mean_freq.loc[target] / mean_freq.loc[background]).sort_values()
    colors = np.where(log2_ratio > 0, BLUE, RED)

    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10, 3.8))
        ax.bar(log2_ratio.index, log2_ratio.to_numpy(), color=colors, width=0.75)
        ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
        ax.set_ylabel(f"log2 ( {target} / {background} )")
        ax.set_title(f"Trinucleotide signature: {target} vs. {background} (canonical k-mers; blue = enriched)")
        ax.tick_params(axis="x", rotation=90, labelsize=8)
        ax.grid(axis="x", visible=False)
        return _save(fig, path)


def plot_read_lengths(reads: pd.DataFrame, palette: Mapping[str, str], path: str | Path) -> Path:
    """Read-length distributions per group (density-normalised step histograms)."""
    bins = np.logspace(np.log10(max(reads["length"].min(), 1)), np.log10(reads["length"].max() + 1), 50)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 3.8))
        for name in _group_order(reads["group"], palette):
            lengths = reads.loc[reads["group"] == name, "length"]
            ax.hist(lengths, bins=bins, density=True, histtype="step", linewidth=2, color=palette[name],
                    label=f"{name} (n={len(lengths):,}, median {lengths.median():,.0f} bp)")
        ax.set_xscale("log")
        ax.set(xlabel="Read length (bp)", ylabel="Density", title="Read length distribution")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(0, -0.16), borderaxespad=0)
        return _save(fig, path)
