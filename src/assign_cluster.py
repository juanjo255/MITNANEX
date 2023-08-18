from psa import psa
from hash_function import hash_table_ids, hash_table_clusters

def assign_cluster(alignment:psa, cluster_pointers:hash_table_ids, clusters:hash_table_clusters) -> int:
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
        new_cluster_id = clusters.set_cluster(alignment)
        cluster_pointers.set_cluster_pointer(alignment.query_id, new_cluster_id)
        cluster_pointers.set_cluster_pointer( alignment.reference_id, new_cluster_id)
