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

    clusters_list = run()

    # I need to plot the coverage of clusters
    clusters_info = pd.DataFrame(
        {
            "coverage": [i.coverage for i in clusters_list.clusters],
            "repr_read_len": [i.longest_read_length for i in clusters_list.clusters],
            "id_longest_read": [i.longest_read_id for i in clusters_list.clusters],
            "id_cluster": [i.id_cluster for i in clusters_list.clusters],
        }
    )

    clusters = clusters_info
    clusters.sort_values(by="id_longest_read", inplace=True)

    convert_fq_to_fa(
        "test/s_cervisae_CEN.PK113-7D_SRR5892449_reads_sample.sorted_reversed.fastq",
        "test/s_cervisae_CEN.PK113-7D_SRR5892449_reads_sample.sorted_reversed.fasta",
    )

    fasta = "test/s_cervisae_CEN.PK113-7D_SRR5892449_reads_sample.sorted_reversed.fasta"
    repr_reads = [i for i in clusters["id_longest_read"]]
    hist = list()
    for i in get_sequences_by_id(fasta, repr_reads):
        ids, seq = i
        hist.append([*count_kmer(k=3, seq=seq).values(), ids])

    ## Get the real mt sequences ##
    with open("test/list_ids_reads_mt.txt", "r") as handle_ids_mt:
        ids_mt = handle_ids_mt.read().splitlines()
        ids_mt = [i[1:] for i in ids_mt]

    hist_df = pd.DataFrame(hist)
    hist_df.rename(columns={hist_df.iloc[:, -1].name: "ids"}, inplace=True)
    hist_df.head()

    ## Dimensionality reduction with PCA and clustering with k-means ##
    pca = PCA(n_components=2)
    pca.fit(hist_df.iloc[:, :-2])

    kmer_reduction = pca.fit_transform(hist_df.iloc[:, :-2])
    kmer_reduction = pd.DataFrame(kmer_reduction)

    ## Merging the dataframe with ids and other relevant information ##
    kmer_reduction = pd.concat([kmer_reduction, hist_df["ids"]], axis=1)
    kmer_reduction.rename(columns={0: "comp1", 1: "comp2"}, inplace=True)
    kmer_reduction = kmer_reduction.merge(
        clusters, how="left", left_on="ids", right_on="id_longest_read"
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
        fasta=fasta, sequences_ids=sequences_ids, output="test/mt_reads_v1.fasta"
    )
