"""Read VCF / depth files with pandas and plot variant-calling results.

Replaces the ``vcfstats`` configuration (``plots_vcf.toml`` + ``custom_macro.py``)
with plain pandas + matplotlib, so no R/plotnine stack is needed.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import plots
from .external import capture
from .io import open_text

log = logging.getLogger(__name__)

TRANSITIONS = {"A>G", "G>A", "C>T", "T>C"}
# Pyrimidine-centred substitution classes (a substitution and its reverse complement are equivalent).
_SUBSTITUTION_CLASS = {
    "C>A": "C>A", "G>T": "C>A", "C>G": "C>G", "G>C": "C>G", "C>T": "C>T", "G>A": "C>T",
    "T>A": "T>A", "A>T": "T>A", "T>C": "T>C", "A>G": "T>C", "T>G": "T>G", "A>C": "T>G",
}


def read_vcf(path: str | Path) -> pd.DataFrame:
    """Load a single-sample VCF into a DataFrame.

    Adds the columns ``DP`` (INFO depth), ``AF`` (allele fraction of the first
    sample's first ALT allele), ``type`` (SNV / insertion / deletion / MNV),
    ``substitution`` (e.g. ``C>T``) and ``pass`` (FILTER == PASS).
    """
    # Count the meta-information lines ('##...') so pandas can start at the header.
    with open_text(path) as handle:
        n_meta = sum(1 for line in handle if line.startswith("##"))
    vcf = pd.read_csv(path, sep="\t", skiprows=n_meta, dtype={"#CHROM": str, "REF": str, "ALT": str})
    vcf = vcf.rename(columns={"#CHROM": "CHROM"})
    if vcf.empty:
        return vcf.assign(DP=[], AF=[], type=[], substitution=[], **{"pass": []})

    # INFO depth: the value of 'DP=' inside the ';'-separated INFO string.
    vcf["DP"] = pd.to_numeric(vcf["INFO"].str.extract(r"(?:^|;)DP=(\d+)")[0], errors="coerce")

    # FORMAT keys + first sample values -> one dict per row -> AF of the first ALT allele.
    sample_col = vcf.columns[9] if len(vcf.columns) > 9 else None
    if sample_col:
        fields = [dict(zip(f.split(":"), v.split(":"))) for f, v in zip(vcf["FORMAT"], vcf[sample_col])]
        vcf["AF"] = pd.to_numeric(pd.Series([f.get("AF", "nan").split(",")[0] for f in fields]), errors="coerce")
    else:
        vcf["AF"] = np.nan

    ref_len, alt = vcf["REF"].str.len(), vcf["ALT"].str.split(",").str[0]
    alt_len = alt.str.len()
    vcf["type"] = np.select(
        [(ref_len == 1) & (alt_len == 1), alt_len > ref_len, alt_len < ref_len], ["SNV", "insertion", "deletion"], "MNV"
    )
    vcf["substitution"] = np.where(vcf["type"] == "SNV", vcf["REF"] + ">" + alt, None)
    vcf["pass"] = vcf["FILTER"] == "PASS"
    return vcf


def read_depth(bam: str | Path) -> pd.DataFrame:
    """Per-position read depth (``samtools depth -a``) as ``CHROM``/``POS``/``depth``."""
    text = capture(["samtools", "depth", "-a", bam])
    return pd.read_csv(io.StringIO(text), sep="\t", names=["CHROM", "POS", "depth"], dtype={"CHROM": str})


def plot_variants_along_genome(vcf: pd.DataFrame, depth: pd.DataFrame, path: str | Path) -> Path:
    """Depth along the mitogenome (top) and allele fraction of each variant (bottom)."""
    palette = {"PASS": plots.BLUE, "filtered": plots.GRAY}
    with plt.rc_context(plots.STYLE):
        fig, (ax_depth, ax_af) = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True,
                                              gridspec_kw={"height_ratios": [1, 1.3]})
        ax_depth.fill_between(depth["POS"], depth["depth"], color=plots.BLUE, alpha=0.25, linewidth=0)
        ax_depth.plot(depth["POS"], depth["depth"], color=plots.BLUE, linewidth=1)
        ax_depth.set(ylabel="Read depth", title=f"Read depth along the mitogenome (median {depth['depth'].median():.0f}x)")

        status = np.where(vcf["pass"], "PASS", "filtered")
        for name in ["filtered", "PASS"]:  # PASS drawn on top
            subset = vcf[status == name]
            ax_af.vlines(subset["POS"], 0, subset["AF"], color=palette[name], linewidth=0.8)
            ax_af.scatter(subset["POS"], subset["AF"], s=22, color=palette[name], edgecolors="white",
                          linewidths=0.5, zorder=3, label=f"{name} (n={len(subset)})")
        ax_af.set_ylim(0, 1.05)
        ax_af.set(xlabel="Position (bp)", ylabel="Allele fraction", title="Variants: allele fraction by position")
        handles, labels = ax_af.get_legend_handles_labels()
        ax_af.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, borderaxespad=0)
        return plots._save(fig, path)


def plot_variant_summary(vcf: pd.DataFrame, path: str | Path) -> Path:
    """Substitution spectrum of PASS SNVs, variant types and allele-fraction distribution."""
    passed = vcf[vcf["pass"]]
    spectrum = (passed["substitution"].dropna().map(_SUBSTITUTION_CLASS).value_counts()
                .reindex(sorted(set(_SUBSTITUTION_CLASS.values())), fill_value=0))
    types = passed["type"].value_counts().reindex(["SNV", "insertion", "deletion", "MNV"], fill_value=0)
    snvs = passed["substitution"].dropna()
    titv = f"Ti/Tv = {snvs.isin(TRANSITIONS).sum()}/{(~snvs.isin(TRANSITIONS)).sum()}" if len(snvs) else ""

    with plt.rc_context(plots.STYLE):
        fig, (ax_spec, ax_type, ax_af) = plt.subplots(1, 3, figsize=(12, 3.6))
        ax_spec.bar(spectrum.index, spectrum.to_numpy(), color=plots.BLUE, width=0.7)
        ax_spec.set(ylabel="PASS SNVs", title=f"Substitution spectrum  {titv}")
        ax_type.bar(types.index, types.to_numpy(), color=plots.BLUE, width=0.7)
        ax_type.set(ylabel="PASS variants", title="Variant types")
        for ax in (ax_spec, ax_type):
            ax.grid(axis="x", visible=False)
        ax_af.hist(passed["AF"].dropna(), bins=np.linspace(0, 1, 21), color=plots.BLUE, edgecolor="white")
        ax_af.set(xlabel="Allele fraction", ylabel="PASS variants",
                  title="Allele fraction (homoplasmy ~1, heteroplasmy <1)")
        return plots._save(fig, path)
