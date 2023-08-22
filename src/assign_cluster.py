from .psa import psa
from .hash_function import hash_table_ids, hash_table_clusters

def assign_cluster(
    alignment: psa, cluster_pointers: hash_table_ids, clusters: hash_table_clusters,
) -> int:
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
        ## If two reads in different cluster, save the read in the cluster of the longest read in the aligment
        ## Otherwise, save both in both clusters
        if query_pointer != reference_pointer:
            

            longest_read_pointer, longest_read_length = (
                (query_pointer, alignment.query_length)
                if alignment.query_length > alignment.reference_length
                else (reference_pointer, alignment.reference_length)
            )
            longest_read_cluster = clusters.get_cluster(longest_read_pointer)
            longest_read_cluster.add_id(alignment.reference_id)
            # Update the longest read of the cluster
            if longest_read_cluster.longest_read < longest_read_length:
                longest_read_cluster.update_longest_read(longest_read_length)
        else:
            # Save in query cluster
            query_cluster = clusters.get_cluster(query_pointer)
            query_cluster.add_id(alignment.reference_id)
            # Update the longest read of the cluster
            if query_cluster.longest_read < alignment.reference_length:
                query_cluster.update_longest_read(alignment.reference_length)

            # Save in reference cluster
            reference_cluster = clusters.get_cluster(reference_pointer)
            reference_cluster.add_id(alignment.query_id)
            # Update the longest read of the cluster
            if reference_cluster.longest_read < alignment.query_length:
                reference_cluster.update_longest_read(alignment.query_length)

    elif query_pointer:
        # Scenario 2: only query_pointer already saved
        query_cluster = clusters.get_cluster(query_pointer)
        query_cluster.add_id(alignment.reference_id)
        # Update the longest read of the cluster
        if query_cluster.longest_read < alignment.reference_length:
            query_cluster.update_longest_read(alignment.reference_length)

    elif reference_pointer:
        # Scenario 3: only query_pointer already saved
        # Add to the existing cluster
        reference_cluster = clusters.get_cluster(reference_pointer)
        reference_cluster.add_id(alignment.query_id)
        # Update the longest read of the cluster
        if reference_cluster.longest_read < alignment.query_length:
            reference_cluster.update_longest_read(alignment.query_length)

    else:
        # Scenario 4: no one saved
        ## We create a new cluster and set the pointer for the reads
        new_cluster_id = clusters.set_cluster(alignment)
        cluster_pointers.set_cluster_pointer(alignment.query_id, new_cluster_id)
        cluster_pointers.set_cluster_pointer(alignment.reference_id, new_cluster_id)

