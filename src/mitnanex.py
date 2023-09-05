from .hash_function import hash_table_ids, hash_table_clusters, estimate_hash_table_size
from .psa import psa
from .assign_cluster import assign_cluster
from .utils import get_sequences_by_id


## PAF FORMAT

# 0	string	query sequence name
# 1	int	query sequence length
# 2	int	query start coordinate (0-based)
# 3	int	query end coordinate (0-based)
# 4	char	‘+’ if query_pointer/target on the same strand; ‘-’ if opposite
# 5	string	Target sequence name
# 6	int	Target sequence length
# 7	int	Target start coordinate on the original strand
# 8	int	Target end coordinate on the original strand
# 9	int	Number of matching bases in the mapping
# 10int	Number bases, including gaps, in the mapping
# 11int	Mapping quality (0-255 with 255 for missing)



# @profile # This is to measure memory consumption
def run() -> hash_table_clusters:
    file = open("test/overlaps_aedes_vexans_all_reads_subsample_sorted_containment.paf", "r")
    #size_table = estimate_hash_table_size('test/aedes_vexans_all_reads_subsample_sorted_length.fastq')
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

    ## Check if a I there are more than one cluster with the same longest read
    # longest_read_per_cluster = [cluster.longest_read_id for cluster in clusters.clusters]
    # unique_longest_reads, count_unique_longest_reads = np.unique(longest_read_per_cluster, return_counts=True)

    # duplicated_clusters = unique_longest_reads[count_unique_longest_reads > 1]
    # print(duplicated_clusters)

    ## Get cluster sequence and cluster
    # repr_reads = [cluster.longest_read_id for cluster in clusters_list.clusters]
    # get_sequences_by_id(repr_reads)
    
    
    return clusters_list


if __name__ == "__main__":
    # root_project_dir = os.path.dirname(os.path.abspath(__name__))
    # sys.path.append(root_project_dir)
    run()
    
