import matplotlib.pyplot as plt
import pandas as pd
from sklearn.neighbors import KernelDensity
import numpy as np
from scipy.signal import argrelextrema
from src.mitnanex import run
import math
import utils_rs
from src.utils import write_fasta
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

prefix = "all_talaro_porechop_18_07_2023"
wd = "test/talaro_santan/mitnanex_results/"

clusters_list = run(f"{wd}{prefix}.paf", 0.7)

# I need to plot the coverage of clusters
coverages_df = pd.DataFrame(
    {
        "coverage": [i.coverage for i in clusters_list.clusters],
        "repr_read_len": [i.longest_read_length for i in clusters_list.clusters],
        "id_longest_read": [i.longest_read_id for i in clusters_list.clusters],
        "id_cluster": [i.id_cluster for i in clusters_list.clusters],
    }
)

clusters = coverages_df.sort_values(by="coverage", ascending=False)
clusters.head()

# # Coverage
#
# The first step is to filter the reads by coverage. For that I will explore the use of Kernel Density since this is clusters_list one dimensional clustering.

clusters["transform"] = [math.log2(i) for i in clusters["coverage"]]
clusters["coverage"].hist(bins=50)
plt.scatter(x=clusters["transform"], y=clusters["coverage"])


kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(
    clusters["coverage"].array.reshape(-1, 1)
)
cov_fdp = kde.score_samples(clusters["coverage"].array.reshape(-1, 1))
plt.plot(clusters["coverage"], cov_fdp)

local_min = argrelextrema(cov_fdp, np.less)[0]
if len(local_min) < 1:
    min_coverage = 3
else:
    min_coverage = clusters.iloc[max(local_min), :]["coverage"]
print("Covertura minima admitida: ", min_coverage)



## ALTERNATIVE TO CALCULATE MINIMUN COVERAGE.
## I WILL LOOK FOR THE GREATES HOLE IN THE DATA, IF THERE IS JUST ONE POINT THEN I CHOOSE THE SECOND GREATEST GAP
def min_cov(clusters):
    cov_gaps = clusters.loc[:, "coverage"].diff(periods=-1).sort_values(ascending=False)
    for k in cov_gaps.index:
        clusters_filt_cov = clusters["coverage"] >= clusters.loc[k, "coverage"]
        if sum(clusters_filt_cov) > 3:
            min_coverage = clusters.loc[k, "coverage"]
            print("Minimum coverage: ", min_coverage)
            return min_coverage

min_coverage = min_cov(clusters)
# FILTER BY COVERAGE
clusters = clusters[clusters["coverage"] >= min_coverage]

## GROUND TRUTH
## Get the real mt sequences
with open("test/talaro_santan/chrMT_reads_ids.txt", "r") as handle_ids_mt:
    ids_mt = handle_ids_mt.read().splitlines()
    ids_mt = [i[:].strip() for i in ids_mt]
clusters["mt"] = [32 if i else 10 for i in clusters["id_longest_read"].isin(ids_mt)]
clusters

# # Oligo composition
# Once I've detected the clusters with higher coverage, 
# which I expect includes mitochondria and contamination, whether nuclear or external, 
# I have to purify these groups. 
# To do this I will use the oligo composition as it is used during metagenomics binning.
# I will create clusters_list script to get the kmers.
# The kmer size will be 3 to solve two things: 
# 1. Intrinsec error from Nanopore
# 2 smaller set as possible (4**3 possible kmers).
# Finally, I will reduce dimensionality using PCA

reads_file = f"{wd}{prefix}_sample.sorted.fastq"
repr_reads = [i for i in clusters["id_longest_read"]]

# RUST VERSION
hist = utils_rs.get_kmer_profiles(repr_reads, reads_file, 3)
hist_df = pd.DataFrame(hist[0])
hist_df["ids"] = hist[1]
hist_df.head()
hist_df.iloc[:, :-1]
hist_df.head()

# # Dimensionality reduction with PCA and clustering with k-means

pca = PCA(n_components=2)
pca.fit(hist_df.iloc[:, :-2])


## Components and variance explained
pca.explained_variance_
# pca.components_


kmer_reduction = pca.fit_transform(hist_df.iloc[:, :-2])
kmer_reduction = pd.DataFrame(kmer_reduction, columns=["comp1", "comp2"])

## merging the dataframe with ids and other relevant information
kmer_reduction["ids"] = hist_df["ids"]
kmer_reduction = kmer_reduction.merge(
    clusters, how="left", left_on="ids", right_on="id_longest_read"
)
kmer_reduction.drop(columns="id_longest_read", inplace=True)

## Annotate which reads are mitochondrial
kmer_reduction["coverage_norm"] = (
    kmer_reduction["coverage"] / kmer_reduction["repr_read_len"]
)
kmer_reduction.head()


sc = plt.scatter(
    kmer_reduction["comp1"],
    kmer_reduction["comp2"],
    c=kmer_reduction["coverage_norm"],
    s=kmer_reduction["mt"],
)
# legend
# Add clusters_list colorbar
cbar = plt.colorbar(sc)
cbar.set_label("Color Scale")
plt.legend()
plt.xlabel("comp1")
plt.ylabel("comp2")


kmeans = KMeans(
    n_clusters=2, max_iter=100, init="k-means++", random_state=0, n_init=1, verbose=1
)
mt_prediction = kmeans.fit_predict(
    kmer_reduction[["comp1", "comp2"]], sample_weight=kmer_reduction["coverage_norm"]
)
kmer_reduction["cluster_prediction"] = mt_prediction


plt.scatter(
    x=kmer_reduction["comp1"],
    y=kmer_reduction["comp2"],
    c=mt_prediction,
    s=kmer_reduction["mt"],
)
plt.xlabel("comp1")
plt.ylabel("comp2")

# # Get the cluster of interest
#
# This step is clusters_list pain in the ass, but since I am looking to keep the free-reference.

selected_cluster_id = (
    kmer_reduction.loc[:, kmer_reduction.columns != "ids"]
    .groupby(by="cluster_prediction")["coverage_norm"]
    .median()
    .idxmax()
)
selected_cluster = kmer_reduction[
    kmer_reduction["cluster_prediction"] == selected_cluster_id
]
selected_cluster.sort_values("coverage_norm", inplace=True, ascending=False)
selected_cluster.head()


sequences_ids = set()
for i in selected_cluster["id_cluster"]:
    sequences_ids.update(clusters_list.get_cluster(i).id_sequences)
write_fasta(
    reads_file=reads_file,
    sequences_ids=sequences_ids,
    output=f"{wd}{prefix}_putative_mt_reads.fasta",
)
