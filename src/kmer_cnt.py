from collections import Counter
from itertools import product
import numpy as np




##### KMER COUNTING #####
## Get complement  complement 
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


## Produce all possible kmers in a sequence using combinations
def possible_kmers(k, ncl='ATCG') -> list:
    kmers = ["".join(i) for i in product(ncl,repeat=k)]
    return kmers


## Get the all kmer from a sequence
## The lexicographically lower kmer between forward and reverse complement is used
## The reason to use the data structure, Counter store hasable elements,
# which make it faster
def count_kmer(k:int, seq:str, id:str):
    ## get the cluster to add the kmer composition
    #cluster = cluster_pointers.get_cluster_pointer(id)
    seq = str(seq)
    len_seq = len(seq)
    counter = Counter(possible_kmers(k=3))
    for i in range(len_seq - k + 1):
        kmer_for = seq[i : (i + k)]
        kmer_rev = kmer_for.translate(comp_tab)[::-1]
        if kmer_for < kmer_rev:
            kmer = kmer_for
        else:
            kmer = kmer_rev
        counter[kmer] += 1
    print(counter)


#count_kmer(3, "AAATTT")