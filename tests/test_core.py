"""Unit tests for the MITNANEX core (no external tools needed): ``pytest tests/``."""

from __future__ import annotations

import gzip
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from mitnanex.clustering import ReadClusterer
from mitnanex.io import read_fasta, read_fastq, read_gfa_segments, read_paf, write_fasta, write_fastq
from mitnanex.kmers import KmerProfiler, gc_content, reverse_complement
from mitnanex.overlaps import OverlapFilter
from mitnanex.selection import CompositionSelector, coverage_threshold_by_gap
from mitnanex.vcf import read_vcf


# --- helpers ------------------------------------------------------------------
def paf_row(q, t, q_len=1000, t_len=1000, qs=0, qe=1000, ts=0, te=1000, strand="+", n_match=900, aln_len=1000):
    return dict(query_id=q, query_len=q_len, query_start=qs, query_end=qe, strand=strand, target_id=t,
                target_len=t_len, target_start=ts, target_end=te, n_match=n_match, aln_len=aln_len, mapq=0)


def naive_canonical_counts(seq: str, k: int) -> Counter:
    """Reference implementation (the logic of the original kmer_cnt.py)."""
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= set("ACGT"):
            counts[min(kmer, reverse_complement(kmer))] += 1
    return counts


# --- io -------------------------------------------------------------------------
def test_fastq_roundtrip_and_id_filter(tmp_path):
    reads = pd.DataFrame({"id": ["r1", "r2", "r3"], "seq": ["ACGT", "GG", "TTTA"], "qual": ["!!!!", '""', "####"]})
    path = write_fastq(reads, tmp_path / "x.fastq")
    loaded = read_fastq(path, with_quality=True)
    assert loaded[["id", "seq", "qual"]].equals(reads)
    assert read_fastq(path, ids=["r3", "r1"])["id"].tolist() == ["r1", "r3"]


def test_fastq_header_with_description_and_gzip(tmp_path):
    path = tmp_path / "x.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@read_1 runid=abc ch=5\nACGT\n+\nIIII\n")
    assert read_fastq(path)["id"].tolist() == ["read_1"]


def test_multiline_fasta(tmp_path):
    seqs = pd.DataFrame({"id": ["a", "b"], "seq": ["ACGT" * 30, "GATTACA"]})
    path = write_fasta(seqs, tmp_path / "x.fa", line_width=50)
    assert read_fasta(path)[["id", "seq"]].equals(seqs)


def test_paf_with_variable_tags(tmp_path):
    path = tmp_path / "x.paf"
    path.write_text("a\t100\t0\t90\t+\tb\t120\t5\t95\t80\t92\t0\ttp:A:S\n"
                    "b\t120\t0\t90\t-\tc\t100\t5\t95\t80\t92\t0\ttp:A:S\tcm:i:5\tdv:f:0.1\n")
    paf = read_paf(path)
    assert paf.shape == (2, 12)
    assert paf["target_id"].tolist() == ["b", "c"]


def test_gfa_segments(tmp_path):
    path = tmp_path / "x.gfa"
    path.write_text("S\tutg1l\tACGTA\tLN:i:5\na\tutg1l\t0\tread1:1-5\t+\t5\nS\tutg2c\tGG\tLN:i:2\n")
    gfa = read_gfa_segments(path)
    assert gfa["id"].tolist() == ["utg1l", "utg2c"] and gfa["length"].tolist() == [5, 2]


# --- overlaps -------------------------------------------------------------------
def test_overlap_classification():
    paf = pd.DataFrame([
        paf_row("a", "b", qs=0, qe=500, ts=500, te=1000),                 # dovetail
        paf_row("a", "c", q_len=500, qs=0, qe=500, ts=200, te=700),       # query contained
        paf_row("a", "d", qs=400, qe=600, ts=400, te=600),                # internal (repeat)
        paf_row("a", "e", qs=0, qe=500, ts=0, te=500, strand="-"),        # reverse-strand dovetail
    ])
    assert OverlapFilter.classify(paf).tolist() == ["dovetail", "query_contained", "internal", "dovetail"]


def test_overlap_filter_identity_and_self():
    paf = pd.DataFrame([paf_row("a", "b", n_match=900), paf_row("a", "c", n_match=300), paf_row("a", "a")])
    kept = OverlapFilter(min_identity=0.6).apply(paf)
    assert kept["target_id"].tolist() == ["b"]


# --- clustering ----------------------------------------------------------------------
def test_greedy_clustering_semantics():
    paf = pd.DataFrame([
        paf_row("r1", "r2", q_len=100, t_len=300),  # new cluster 0 {r1, r2}
        paf_row("r3", "r2", q_len=500, t_len=300),  # r3 joins cluster 0
        paf_row("r4", "r5", q_len=50, t_len=60),    # new cluster 1 {r4, r5}
        paf_row("r1", "r4", q_len=100, t_len=50),   # both clustered -> no merge
    ])
    result = ReadClusterer().fit(paf)
    clusters = result.clusters.set_index("cluster_id")
    assert clusters.loc[0, "coverage"] == 3 and clusters.loc[0, "repr_read_id"] == "r3"
    assert clusters.loc[1, "coverage"] == 2 and clusters.loc[1, "repr_read_id"] == "r5"
    assert set(result.reads_of([1])) == {"r4", "r5"}


# --- k-mers ------------------------------------------------------------------------------
@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_kmer_counts_match_naive(k):
    rng = np.random.default_rng(k)
    seq = "".join(rng.choice(list("ACGTN"), size=500, p=[0.24, 0.24, 0.24, 0.24, 0.04]))
    profiler = KmerProfiler(k)
    expected = naive_canonical_counts(seq, k)
    counts = dict(zip(profiler.canonical_kmers, profiler.counts(seq)))
    assert counts == {kmer: expected.get(kmer, 0) for kmer in profiler.canonical_kmers}


def test_pooled_profiles_weight_reads_by_length():
    profiler = KmerProfiler(1)
    seqs = pd.Series(["AAAA", "CC"])  # pooled: 4 A + 2 C over 6 bases
    pooled = profiler.pooled_profiles(seqs, pd.Series(["x", "x"]))
    assert pooled.loc["x", "A"] == pytest.approx(4 / 6) and pooled.loc["x", "C"] == pytest.approx(2 / 6)


def test_kmer_profile_strand_independent():
    profiler = KmerProfiler(3)
    seq = "ACGTTGCAAGGCTTAACCGGTA" * 10
    assert len(profiler.canonical_kmers) == 32
    np.testing.assert_allclose(profiler.profile(seq), profiler.profile(reverse_complement(seq)))
    assert profiler.profile(seq).sum() == pytest.approx(1)


def test_gc_content():
    assert gc_content(pd.Series(["GGCC", "ATAT", "acgt"])).tolist() == [1.0, 0.0, 0.5]


# --- selection ------------------------------------------------------------------------------
def test_coverage_gap_threshold():
    coverage = pd.Series([500, 480, 470, 460, 60, 50, 45, 40, 3])
    assert coverage_threshold_by_gap(coverage, min_clusters=3) == 460
    # Not enough clusters above any gap -> keep everything.
    assert coverage_threshold_by_gap(pd.Series([10, 5]), min_clusters=3) == 5


def _composition_case(rng):
    """Toy candidates: one high-coverage cluster with mt-like composition (AT-rich),
    one GC-rich repeat cluster and three low-coverage nuclear clusters."""
    def random_seq(gc, n=3000):
        p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        return "".join(rng.choice(list("ACGT"), size=n, p=p))
    seqs = {"mt": random_seq(0.44), "rdna": random_seq(0.58), "n1": random_seq(0.40),
            "n2": random_seq(0.40), "n3": random_seq(0.41)}
    candidates = pd.DataFrame({
        "cluster_id": range(5), "coverage": [500, 230, 57, 55, 50],
        "repr_read_id": list(seqs), "repr_read_len": [16000, 30000, 24000, 21000, 32000],
    })
    cluster_seqs = pd.DataFrame({"cluster_id": range(5), "seq": list(seqs.values())})
    return candidates, cluster_seqs


def test_group_rule_max_keeps_highest_coverage_cluster():
    candidates, seqs = _composition_case(np.random.default_rng(0))
    table = CompositionSelector(group_rule="max").fit(candidates, seqs).table
    assert table.set_index("repr_read_id").loc["mt", "selected"]


def test_group_rule_median_reproduces_original_failure():
    candidates, seqs = _composition_case(np.random.default_rng(0))
    table = CompositionSelector(group_rule="median").fit(candidates, seqs).table.set_index("repr_read_id")
    # The original rule picks the GC-distinct repeat over the mitochondrial cluster.
    assert table.loc["rdna", "selected"] and not table.loc["mt", "selected"]


# --- vcf -------------------------------------------------------------------------------------------
def test_read_vcf(tmp_path):
    path = tmp_path / "x.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "chrM\t73\t.\tA\tG\t.\tPASS\tDP=100;TLOD=50\tGT:AD:AF\t0/1:2,98:0.98\n"
        "chrM\t310\t.\tT\tTC\t.\tstrand_bias\tDP=80\tGT:AD:AF\t0/1:60,20:0.25\n"
    )
    vcf = read_vcf(path)
    assert vcf["type"].tolist() == ["SNV", "insertion"]
    assert vcf["AF"].tolist() == [0.98, 0.25] and vcf["DP"].tolist() == [100, 80]
    assert vcf["pass"].tolist() == [True, False] and vcf.loc[0, "substitution"] == "A>G"


# --- assembly (needs minimap2) --------------------------------------------------------------------
@pytest.mark.skipif(__import__("shutil").which("minimap2") is None, reason="minimap2 not installed")
def test_collapse_circular_dimer():
    from mitnanex.assembly import collapse_tandem_copies

    rng = np.random.default_rng(3)
    unit = "".join(rng.choice(list("ACGT"), size=20_000))
    rotated_dimer = unit[5_000:] + unit + unit[:5_000]  # Flye dimers start anywhere in the circle
    contigs = pd.DataFrame({"id": ["contig_1"], "seq": [rotated_dimer], "length": [len(rotated_dimer)]})
    collapsed, period = collapse_tandem_copies(contigs)
    assert period == len(unit) and collapsed.loc[0, "length"] == len(unit)

    monomer = pd.DataFrame({"id": ["contig_1"], "seq": [unit], "length": [len(unit)]})
    assert collapse_tandem_copies(monomer)[1] is None
