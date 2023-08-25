from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt


## I need to separate clusters of reads by coverage
## NOTE: so far only no plants organisms will be tackle.


file = pd.read_csv("overlaps_talaro_18_07_2023.paf",delimiter='\t',header=None)
file [10] = file[10].astype(int)
file = file[file[10]>500]
file = (file.groupby(by=0))
## This are the ids to get the reads with seqkit
#clusters = pd.DataFrame(file.size().index)
#clusters.to_csv('ids.txt',sep='\t',index=False)

reads_per_cluster_gc = pd.read_csv('gc_content_per_cluster.tsv', delimiter='\t',header=None)
reads_per_cluster_gc.sort_values(by=0,inplace=True)
reads_per_cluster_gc['covertura2'] = 0#file.size().values
reads_per_cluster_gc['covertura'] = file.size().values
reads_per_cluster_gc.rename(columns={2:'gc_content'}, inplace=True)

## First view before Kmeans
#plt.scatter(x=reads_per_cluster_gc['gc_content'], y=reads_per_cluster_gc['covertura'])
#plt.show()

kmeans = KMeans(n_clusters=20, max_iter=100, init='k-means++', random_state=0, n_init=1, verbose=1)
prediction = kmeans.fit_predict(reads_per_cluster_gc.loc[:, ['covertura','gc_content']])


plt.scatter(x=reads_per_cluster_gc['gc_content'], y=reads_per_cluster_gc['covertura'], c=prediction)

plt.show()
