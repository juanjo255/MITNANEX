
# Add the root project directory to the Python path
import sys
import os
root_project_dir = os.path.dirname(os.path.abspath(__name__))
sys.path.append(root_project_dir)
from src.hash_function import hash_table_ids, hash_table_clusters
from src.psa import psa
from src.assign_cluster import assign_cluster




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



def filter_alignment(alignment: str) -> psa:
    """Take one line of a paf format given by minimap2, filter it and save it in a cluster
    Args:
        align (str): One line from paf file which contains the default paf format given
            by ava-ont function in minimap2.

    Returns:
        cluster: class where the read will be gathered otherwise it will initialized a new cluster.
    """
    alignment = psa(alignment.split("\t"))
    ## filter mappings shorter than 500 pb
    length_threshold = 500
    if (int(alignment.align_length) < length_threshold):
        return None

    ## filter only hits that where one read contains the other one equal or greater to threshold
    ## NOTE: So far I do not take into account the difference between internal matches and containments
    ## Here containment is just the read aligned to other read at least threshold_contaiment
    threshold_containment = 0.8
    if alignment.align_length >= min(alignment.query_length, alignment.reference_length) * threshold_containment and (alignment.match_bases/alignment.align_length) > 0.6:
        return alignment

    return None 


# @profile # This is to measure memory consumption
def run () -> hash_table_clusters:
    file = open("test/aedes_vexans_reads_subsample_overlaps_filtered.paf", "r")
    alignment = file.readline().strip()
    clusters = hash_table_clusters()
    cluster_pointers = hash_table_ids(size_table=int(1e10))
    
    while alignment:
    
        #print(alignment)
        #print("Filtrando alignment")
        pass_alignment = filter_alignment(alignment)
        if pass_alignment:
            #print("Asignando cluster")
            assign_cluster(pass_alignment, cluster_pointers, clusters)

        # New align
        alignment = file.readline().strip()

    file.close()

    ## Check if a I there are more than one cluster with the same longest read
    # longest_read_per_cluster = [cluster.longest_read_id for cluster in clusters.clusters]
    # unique_longest_reads, count_unique_longest_reads = np.unique(longest_read_per_cluster, return_counts=True) 
    
    # duplicated_clusters = unique_longest_reads[count_unique_longest_reads > 1]
    # print(duplicated_clusters)
    return clusters
    
if __name__ == "__main__":
    run()
    
