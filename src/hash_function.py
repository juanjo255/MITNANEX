import numpy as np
from .cluster import cluster
from .psa import psa
from collections import Counter
## This function is to do something like a hash table
## where each id is convert it to a position in an array

## NOTE: This arrays does not contain the cluster of reads,
## it contains a pointer to other array which contain the objects where reads are contain
class hash_table_ids:
    
    def __init__(self, size_table=int(1e10)) -> None:
        # self.size_table = size_table
        # self.read_ids = np.zeros(size_table, dtype=np.uint32)
        self.read_ids = Counter()

    def set_cluster_pointer (self, key:str, value:int) -> None:
        #hash_key = hash(key) % self.size_table
        self.read_ids[key] = value

    def get_cluster_pointer (self, key:str) -> int:
        #hash_key = hash(key) % self.size_table
        
        return self.read_ids[key]

class hash_table_clusters:
    
    def __init__(self) -> None:
        self.clusters = []

    def set_cluster (self, alignment:psa) -> int:
        new_cluster = cluster(alignment)
        self.clusters.append(new_cluster)
        
        ## NOTE: This is not the correct index but is to avoid the zero. 
        ## Since hash_table_ids initialize an array of zeros.
        ## The idea is to link the index of the clusters in clusters with the reads in hash_table_ids
        new_cluster.set_cluster_id(len(self.clusters))
        return len(self.clusters)

    def get_cluster (self, key:int) -> cluster:
        ## NOTE: Check set_cluster method
        return self.clusters[key-1]

def estimate_hash_table_size(path:str,path2) -> int:
    """number of reads in fast file

    Args:
        path (str): location of the file
    """
    file = open(path, "r")
    file2 = open(path2, "r")
    line = file.readline().strip()
    line2 = file2.readline().strip()
    size = 0
    size2 = 0
    while line:
        # if line.startswith(('A','T','G','C')):
        #     size+=1
        if line.split("\t")[0] != file.readline().strip().split("\t")[0]:
            #print(line.split("\t")[0])
            size+=1
        # New alignment
        line = file.readline().strip()
    while line2:
        if line2.startswith(('A','T','G','C')):
            size2+=1
        # New alignment
        line2 = file2.readline().strip()
    file.close()
    file2.close()
    print(size, size2)
    return size


