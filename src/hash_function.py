import numpy as np

## This function is to do something like a hash table
## where each id is convert it to a position in an array

## NOTE: This arrays does not contain the cluster of reads,
## it contains a pointer to other array which contain the objects where reads are contain
class hash_ids():
    
    def __init__(self, size_table=int(1e5)) -> None:
        self.size_table = np.zeros(size_table, dtype=int)

    def __setitem__ (self, hash_key, value):
        self.size_table[hash_key] = value

    def __getitem__(self, hash_key):
        return self.size_table[hash_key]

a = hash_ids()

