from src.mitnanex import run
from src.utils import write_fasta
from src.kmer_reduction_PCA import kmer_reduction
from src.cluster_kmer_profiles import cluster_kmer_profiles
from src.minimun_cov import set_minimun_cov
import utils_rs
import pandas as pd
import sys


if __name__ == "__main__":
    ## GET INPUTS
    args = sys.argv
    reads_file = args[1]  # sample read file
    paf_file = args[2]  # PAF
    coverage = int(args[3])  # minimun coverage per cluster
    map_identity = float(args[4])  # minimun coverage per cluster
    output = args[5]  # file to write mt reads

    # MAIN PROGRAM
    clusters_list = run(paf_file, map_identity)

    ## Gather clusters info
    clusters_info = pd.DataFrame(
        {
            "coverage": [i.coverage for i in clusters_list.clusters],
            "repr_read_len": [i.longest_read_length for i in clusters_list.clusters],
            "id_longest_read": [i.longest_read_id for i in clusters_list.clusters],
            "id_cluster": [i.id_cluster for i in clusters_list.clusters],
        }
    )
    
    ## Get minimum coverage
    min_coverage = set_minimun_cov(clusters_info, coverage)

    ## Filter clusters by coverage
    clusters_info = clusters_info[clusters_info["coverage"] >= min_coverage]
    clusters_info["coverage_norm"] = (
        clusters_info["coverage"] / clusters_info["repr_read_len"]
    )

    ## Get kmer composition from representative reads from all the cluster passed
    repr_reads = [i for i in clusters_info["id_longest_read"]]
    kmer_profiles, ids = utils_rs.get_kmer_profiles(repr_reads, reads_file, 3)
    kmer_profiles_df = pd.DataFrame(kmer_profiles)
    kmer_profiles_df["ids"] = ids

    ## Dimensionality reduction with PCA ##
    ## Additionally we merge information from cluster_info
    kmer_reduction_df = kmer_reduction(
        kmer_profiles_df=kmer_profiles_df, clusters_info=clusters_info, n_comp=2
    )

    ## Clustering reads using K-means expecting 2 clusters ##
    kmer_reduction_df["cluster_prediction"] = cluster_kmer_profiles(
        kmer_reduction_df=kmer_reduction_df, n_clusters=2, max_iter=100
    )

    ## Get the cluster of interest ##
    # This step is a pain in the ass,
    # but since I am looking to keep the free-reference,
    # I will select the cluster with the highest average coverage.
    selected_cluster_id = (
        kmer_reduction_df.loc[:, kmer_reduction_df.columns != "ids"]
        .groupby(by="cluster_prediction")["coverage_norm"]
        .median()
        .idxmax()
    )
    selected_cluster = kmer_reduction_df[
        kmer_reduction_df["cluster_prediction"] == selected_cluster_id
    ]
    # print(selected_cluster)
    # import matplotlib.pyplot as plt
    # plt.scatter(x=kmer_reduction_df['comp1'], y=kmer_reduction_df['comp2'], c=kmer_reduction_df['cluster_prediction'])
    # plt.xlabel('comp1')
    # plt.ylabel('comp2')
    # plt.show()

    ## Get sequences from selected clusters and write fasta
    sequences_ids = set()
    for i in selected_cluster["id_cluster"]:
        sequences_ids.update(clusters_list.get_cluster(i).id_sequences)

    write_fasta(
        reads_file=reads_file,
        sequences_ids=sequences_ids,
        output=output,
    )
