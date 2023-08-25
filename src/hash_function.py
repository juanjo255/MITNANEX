import numpy as np
from src.cluster import cluster
from src.psa import psa
## This function is to do something like a hash table
## where each id is convert it to a position in an array

## NOTE: This arrays does not contain the cluster of reads,
## it contains a pointer to other array which contain the objects where reads are contain
class hash_table_ids:
    
    def __init__(self, size_table=int(1e10)) -> None:
        self.size_table = size_table
        self.read_ids = np.zeros(size_table, dtype=np.uint32)

    def set_cluster_pointer (self, key:str, value:int) -> None:
        hash_key = hash(key) % self.size_table
        self.read_ids[hash_key] = value

    def get_cluster_pointer (self, key:str) -> int:
        hash_key = hash(key) % self.size_table
        
        return self.read_ids[hash_key]

class hash_table_clusters:
    
    def __init__(self) -> None:
        self.clusters = []

    def set_cluster (self, alignment:psa) -> int:
        new_cluster = cluster(alignment)
        self.clusters.append(new_cluster)
        
        ## NOTE: This is not the correct index but is to avoid the zero. 
        ## Since hash_table_ids initialize an array of zeros.
        ## The idea is to link the index of the clusters in clusters with the reads in hash_table_ids
        
        return len(self.clusters)

    def get_cluster (self, key:int) -> cluster:
        ## NOTE: Check set_cluster method
        return self.clusters[key-1]

