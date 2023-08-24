from sklearn.cluster import KMeans
import pandas as pd


## I need to separate clusters of reads by coverage
## NOTE: so far only no plants organisms will be tackle.


file = pd.read_csv("overlaps_talaro_18_07_2023.paf",delimiter='\t',header=None)
file [10] = file[10].astype(int)
file = file[file[10]>500]
file = (file.groupby(by=0))#.size()
print(file.size().sort_values(ascending=False).describe())
#unique_readsfile.size().index


# gc_content = pd.read_csv('../longest_read_every_cluster_gc_content.tsv', delimiter='\t',header=None)
# kmeas = KMeans(n_clusters=2, max_iter=100, init='k-means++', random_state=0).fit(gc_content[2])
#print(gc_content[2].head())
