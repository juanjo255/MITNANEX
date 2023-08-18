from hash_function import hash_table_ids, hash_table_clusters
from psa import psa


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


def asign_cluster(alignment:psa, cluster_pointers:hash_table_ids, clusters:hash_table_clusters) -> int:
    """Asign a cluster to a pair of reads

    Args:
        alignment_fields (list): alignment_fields of a pairwise alignment in PAF format
        cluster_pointers (hash_table_ids): Contains pointers to clusters
        clusters (hash_table_clusters): Contains clusters

    Returns:
        int: 
    """
    
    ## Get the cluster pointer if it exists in table for both reads query_pointer and reference_pointer
    query_pointer = cluster_pointers.get_cluster_pointer(alignment.query_id)
    reference_pointer = cluster_pointers.get_cluster_pointer(alignment.reference_id)

    if query_pointer and reference_pointer:
        # Scenario 1: both of them already saved
        # NOTE: Can two reads in two different clusters? 
        pass
    elif query_pointer:
        # Scenario 2: only query_pointer already saved
        query_cluster = clusters.get_cluster(query_pointer)
        query_cluster.add_id(alignment.reference_id)
    elif reference_pointer:
        # Scenario 3: only query_pointer already saved
        # Add to the existing cluster
        reference_cluster = clusters.get_cluster(reference_pointer)
        reference_cluster.add_id(alignment.query_id)
    else:
        # Scenario 4: no one saved
        ## We create a new cluster and set the pointer for the reads
        new_cluster = clusters.set_cluster()
        cluster_pointers.set_cluster_pointer(new_cluster, alignment.query_id)
        cluster_pointers.set_cluster_pointer(new_cluster, alignment.reference_id)

def filter_alignment(alignment: psa, cluster_pointers:hash_table_ids):
    """Take one line of a paf format given by minimap2, filter it and save it in a cluster
    Args:
        align (str): One line from paf file which contains the default paf format given
            by ava-ont function in minimap2.

    Returns:
        cluster: class where the read will be gathered otherwise it will initialized a new cluster.
    """

    ## filter mappings shorter than 500 pb
    length_threshold = 500
    if (int(alignment.align_length) < length_threshold):
        return None

    ## filter only hits that where one read contains the other one equal or greater to threshold
    ## NOTE: So far I do not take into account the difference between internal matches and containments
    ## Here containment is just the read aligned to other read at least threshold_contaiment
    threshold_containment = 0.8
    if alignment.align_length >= min(alignment.query_length, alignment.refence_length) * threshold_containment:
        pass


# @profile # This is to measure memory consumption
def run () -> None:
    file = open("overlaps_talaro_18_07_2023.paf", "r")
    alignment = file.readline()
    # while alignment:
    #    alignment = file.readline()
    
    cluster_pointers = hash_table_ids(size_table=int(1e10))
    clusters = hash_table_clusters
    alignment = psa(alignment.split("\t"))
    filter_alignment(alignment, cluster_pointers, clusters)
    
    
if __name__ == "__main__":
    chunksize = 1000
    run()
    
