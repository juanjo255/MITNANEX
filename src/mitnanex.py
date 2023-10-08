from .hash_function import hash_table_ids, hash_table_clusters
from .psa import psa
from .assign_cluster import assign_cluster
from .utils import write_fasta
from .kmer_reduction_PCA import kmer_reduction
from .cluster_kmer_profiles import cluster_kmer_profiles
from .minimun_cov import set_minimun_cov
import utils
import pandas as pd
import sys
import math

## PAF FORMAT

# 0	string	query sequence name
# 1	int	query sequence length
# 2	int	query start coordinate (0-based)
# 3	int	query end coordinate (0-based)
# 4	char ‘+’ if query_pointer/target on the same strand; ‘-’ if opposite
# 5	string	Target sequence name
# 6	int	Target sequence length
# 7	int	Target start coordinate on the original strand
# 8	int	Target end coordinate on the original strand
# 9	int	Number of matching bases in the mapping
# 10int	Number bases, including gaps, in the mapping
# 11int	Mapping quality (0-255 with 255 for missing)


# @profile # This is to measure memory consumption
def run(paf: str) -> hash_table_clusters:
    file = open(paf, "r")
    alignment = file.readline().strip()
    clusters_list = hash_table_clusters()
    cluster_pointers = hash_table_ids(size_table=int(1e5))

    # Iterate through each alignment
    while alignment:
        alignment = psa(alignment.strip().split("\t"))
        if alignment.map_identity >= 0.6:
            assign_cluster(alignment, cluster_pointers, clusters_list)
        # New alignment
        alignment = file.readline().strip()

    file.close()
    return clusters_list


if __name__ == "__main__":
    ## GET INPUTS
    args = sys.argv
    reads_file = args[1]  # sample read file
    paf_file = args[2]  # PAF
    coverage = args[3]  # minimun coverage per cluster
    output = args[4]  # file to write mt reads

    # MAIN PROGRAM
    clusters_list = run(
        paf_file
    )  # "test/sara_reads/jfdminion11_sara_reads_containments.paf"

    ## Gather clusters info
    clusters_info = pd.DataFrame(
        {
            "coverage": [i.coverage for i in clusters_list.clusters],
            "repr_read_len": [i.longest_read_length for i in clusters_list.clusters],
            "id_longest_read": [i.longest_read_id for i in clusters_list.clusters],
            "id_cluster": [i.id_cluster for i in clusters_list.clusters],
        }
    )

    ## Normalize coverage

    ## Get minimum coverage
    min_coverage = set_minimun_cov(clusters_info, coverage)

    ## Filter clusters by coverage
    clusters_info = clusters_info[clusters_info['coverage'] >= min_coverage]
    clusters_info['coverage_norm'] = clusters_info['coverage'] / clusters_info['repr_read_len']
    
    ## Get kmer composition from representative reads from all the cluster passed
    repr_reads = [i for i in clusters_info["id_longest_read"]]
    kmer_profiles, ids = utils.get_kmer_profiles(repr_reads, reads_file, 3)
    kmer_profiles_df = pd.DataFrame(kmer_profiles)
    kmer_profiles_df["ids"] = ids

    ## Dimensionality reduction with PCA ##
    ## Additionally we merge information from cluster_info
    kmer_reduction_df = kmer_reduction(
        kmer_profiles_df=kmer_profiles_df, clusters_info=clusters_info, n_comp=2
    )

    ## Clustering reads using K-means expecting 2 clusters ##
    kmer_reduction_df["cluster_prediction"] = cluster_kmer_profiles(
        kmer_reduction_df=kmer_profiles_df, n_clusters=2, max_iter=100
    )

    ## Get the cluster of interest ##
    # This step is a pain in the ass,
    # but since I am looking to keep the free-reference,
    # I will select the cluster with the highest average coverage.
    selected_cluster_id = (
        kmer_reduction.loc[:, kmer_reduction.columns != "ids"]
        .groupby(by="cluster_prediction")["coverage"]
        .median()
        .idxmax()
    )
    selected_cluster = kmer_reduction[
        kmer_reduction["cluster_prediction"] == selected_cluster_id
    ]

    ## Get sequences from selected clusters and write fasta
    sequences_ids = set()
    for i in selected_cluster["id_cluster"]:
        sequences_ids.update(clusters_list.get_cluster(i).id_sequences)

    write_fasta(
        reads_file=reads_file,
        sequences_ids=sequences_ids,
        output=output,
    )
