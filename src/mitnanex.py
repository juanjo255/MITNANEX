from .hash_function import hash_table_ids, hash_table_clusters
from .psa import psa
from .assign_cluster import assign_cluster
from .utils import convert_fq_to_fa, write_fasta
from .kmer_cnt import get_kmer_profiles
from .kmer_reduction_PCA import kmer_reduction
from .cluster_kmer_profiles import cluster_kmer_profiles
import pandas as pd

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
def run(paf:str) -> hash_table_clusters:
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
    # root_project_dir = os.path.dirname(os.path.abspath(__name__))
    # sys.path.append(root_project_dir)

    # FIXME: Is this true?
    # Transform fastq to reads_file so it is faster to iterate 
    reads_file = "test/sara_reads/jfdminion11_sara_reads_sample.sorted.fastq"

    # MAIN PROGRAM
    clusters_list = run("test/sara_reads/jfdminion11_sara_reads_containments.paf")

    ## Gather clusters info
    clusters_info = pd.DataFrame(
        {
            "coverage": [i.coverage for i in clusters_list.clusters],
            "repr_read_len": [i.longest_read_length for i in clusters_list.clusters],
            "id_longest_read": [i.longest_read_id for i in clusters_list.clusters],
            "id_cluster": [i.id_cluster for i in clusters_list.clusters],
        }
    )
    clusters_info.sort_values(by="id_longest_read", inplace=True)

    ## Get kmers from representative reads from each cluster
    repr_reads = [i for i in clusters_info["id_longest_read"]]
    kmer_profiles_df = get_kmer_profiles (repr_reads, reads_file)

    ## Dimensionality reduction with PCA and clustering with k-means ##
    kmer_reduction_df = kmer_reduction(kmer_profiles_df, clusters_info)

    ## Clustering reads using K-means expecting 2 clusters ##
    cluster_prediction = cluster_kmer_profiles(kmer_profiles_df)
    kmer_reduction["cluster_prediction"] = cluster_prediction

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
        reads_file=reads_file, sequences_ids=sequences_ids, output="test/sara_reads/mt_reads_v1.fasta"
    )
