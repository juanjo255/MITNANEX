from .hash_function import hash_table_ids, hash_table_clusters
from .psa import psa
from .assign_cluster import assign_cluster

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
def run(paf: str, map_identity:float) -> hash_table_clusters:
    file = open(paf, "r")
    alignment = file.readline().strip()
    clusters_list = hash_table_clusters()
    cluster_pointers = hash_table_ids(size_table=int(1e5))

    # Iterate through each alignment
    while alignment:
        alignment = psa(alignment.strip().split("\t"))
        if alignment.map_identity >= map_identity:
            ## filter only hits that where one read contains the other one equal or greater to threshold
            ## NOTE: So far I do not take into account the difference between internal matches and containments
            ## Here containment is just the read aligned to other read at least threshold_contaiment
            threshold_containment = 0.8
            if alignment.align_length >= min(alignment.query_length, alignment.reference_length) * threshold_containment:
                assign_cluster(alignment, cluster_pointers, clusters_list)
                
        # New alignment
        alignment = file.readline().strip()

    

    file.close()
    return clusters_list