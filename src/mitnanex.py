from .hash_function import hash_table_ids, hash_table_clusters
from .psa import psa
from .assign_cluster import assign_cluster
from .utils import get_sequences_by_id, convert_fq_to_fa, write_fasta
from .kmer_cnt import count_kmer
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

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
def run() -> hash_table_clusters:
    file = open("test/s_cervisae_CEN.PK113-7D_containments_sorted_reversed.paf", "r")
    alignment = file.readline().strip()
    clusters_list = hash_table_clusters()
    cluster_pointers = hash_table_ids(size_table=int(1e5))

    # Iterate through each alignment
    while alignment:
        alignment = psa(alignment.strip().split("\t"))
        if alignment.map_identity >= 0.5:
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
    reads_file = "test/s_cervisae_CEN.PK113-7D_SRR5892449_reads_sample.sorted_reversed.fastq"
    reads_file= convert_fq_to_fa(
        fastq=reads_file,
        output= "".join(reads_file.split(".")) + '.fasta',
    )

    # MAIN PROGRAM
    clusters_list = run()

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
    kmer_profiles = list()
    for i in get_sequences_by_id(reads_file, repr_reads):
        ids, seq = i
        kmer_profiles.append([*count_kmer(k=3, seq=seq).values(), ids])
    kmer_profiles_df = pd.DataFrame(kmer_profiles)
    kmer_profiles_df.rename(columns={kmer_profiles_df.iloc[:, -1].name: "ids"}, inplace=True)

    ## Dimensionality reduction with PCA and clustering with k-means ##
    pca = PCA(n_components=2)
    kmer_reduction = pd.DataFrame(pca.fit_transform(kmer_profiles_df.iloc[:, :-2]), columns=['comp1','comp2'])

    ## Merging the dataframe with ids and other relevant information ##
    kmer_reduction['ids'] = kmer_profiles_df["ids"]
    kmer_reduction = kmer_reduction.merge(
        clusters_info, how="left", left_on="ids", right_on="id_longest_read"
    )
    kmer_reduction.drop(columns="id_longest_read", inplace=True)

    ## Filter reads by coverage ##
    kmer_reduction = kmer_reduction[kmer_reduction["coverage"] > 10]

    ## Clustering reads using K-means expecting 2 clusters ##
    kmeans = KMeans(
        n_clusters=2,
        max_iter=100,
        init="k-means++",
        random_state=0,
        n_init=1,
        verbose=1,
    )
    mt_prediction = kmeans.fit_predict(kmer_reduction[["comp1", "comp2"]])
    kmer_reduction["cluster_prediction"] = mt_prediction

    ## Get the cluster of interest ##
    # This step is a pain in the ass,
    # but since I am looking to keep the free-reference,
    # I will select the cluster with the highest average coverage.
    selected_cluster_id = (
        kmer_reduction.loc[:, kmer_reduction.columns != "ids"]
        .groupby(by="cluster_prediction")["coverage"]
        .mean()
        .idxmax()
    )
    selected_cluster = kmer_reduction[
        kmer_reduction["cluster_prediction"] == selected_cluster_id
    ]

    sequences_ids = set()
    for i in selected_cluster["id_cluster"]:
        sequences_ids.update(clusters_list.get_cluster(i).id_sequences)
        
    write_fasta(
        reads_file=reads_file, sequences_ids=sequences_ids, output="test/mt_reads_v1.reads_file"
    )
