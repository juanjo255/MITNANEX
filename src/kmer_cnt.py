from collections import Counter
from itertools import product
import pandas as pd
from .utils import get_sequences_by_id

##### KMER COUNTING #####
## Get reverse complement
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


## Produce all possible kmers in a sequence using combinations
def possible_kmers(k, ncl="ATCG") -> list:
    kmers = ["".join(i) for i in product(ncl, repeat=k)]
    canon_kmers = [
        kmer for kmer in kmers if kmer == min(kmer, kmer.translate(comp_tab)[::-1])
    ]
    return canon_kmers


## Get the all kmer from a sequence
## The lexicographically lower kmer between forward and reverse complement is used
## The reason to use the data structure, Counter store hasable elements,
# which make it faster
def count_kmer(k: int, seq: str ) -> Counter:
    ## get the cluster to add the kmer composition
    # cluster = cluster_pointers.get_cluster_pointer(id)
    possible_kmers_list: list=possible_kmers(k)
    seq = str(seq)
    len_seq = len(seq)
    tot_kmer_expected = len_seq - k + 1
    kmer_counter = Counter(possible_kmers_list)
    for i in range(tot_kmer_expected):
        kmer_for = seq[i : (i + k)]
        kmer_rev = kmer_for.translate(comp_tab)[::-1]
        if kmer_for < kmer_rev:
            kmer = kmer_for
        else:
            kmer = kmer_rev

        kmer_counter[kmer] += 1
    kmer_counter = norm_kmers(kmer_counter, tot_kmer_expected)
    return kmer_counter

def norm_kmers (kmer_counter:Counter, tot_kmer_expected:int) -> Counter:
    for kmer in kmer_counter:
        if kmer_counter[kmer] > 1:
            kmer_counter[kmer] = (kmer_counter[kmer] - 1)  / tot_kmer_expected
        else:
            kmer_counter[kmer] = float(kmer_counter[kmer])
    return kmer_counter

# possible_kmers_list = possible_kmers(k=3)
# print(len(possible_kmers_list))
# count_kmer(3, "AAATTT", possible_kmers_list)
# print(len(possible_kmers_list))

def get_kmer_profiles (repr_reads:list, reads_file:str) -> pd.DataFrame:
        kmer_profiles = list()
        for i in get_sequences_by_id(reads_file, repr_reads):
            ids, seq = i
            kmer_profiles.append([*count_kmer(k=3, seq=seq).values(), ids])
        kmer_profiles_df = pd.DataFrame(kmer_profiles)
        kmer_profiles_df.rename(columns={kmer_profiles_df.iloc[:, -1].name: "ids"}, inplace=True)
        return kmer_profiles_df