import time
from src.mitnanex import run
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans

inicio = time.time()
a = run()

# I need to plot the coverage of clusters 
coverages_df = pd.DataFrame({"coverage":[i.coverage for i in a.clusters], "longest_read_len":[i.longest_read_length for i in a.clusters], "id_longest_read": [i.longest_read_id for i in a.clusters] })
#coverages_df.sort_values(by='coverage',inplace=True, ascending=False)
coverages_df ['GC_percentage'] = 0
file = open('test/aedes_vexans_mt_reads_subsample.fastq', 'r')
oneline = file.readline()
while oneline:
    if oneline.startswith("@") and len(oneline.split(" ")) > 1:
        id = (oneline.split(" ")[0].strip())[1:]
        if id in coverages_df['id_longest_read'].values:
            # get ADN (it's the line following the id line)
            #print(id)
            sequence=file.readline()
            cg_count = sequence.count('C') + sequence.count("G")
            coverages_df.loc[coverages_df["id_longest_read"]==id,'GC_percentage'] = cg_count/len(sequence)
            
    oneline=file.readline()

print(coverages_df.describe())
final = time.time()
tiempo = final - inicio
print("my code", tiempo)

## First view before Kmeans
#plt.scatter(x=coverages_df['GC_percentage'], y=coverages_df['coverage'])

kmeans = KMeans(n_clusters=3, max_iter=100, init='k-means++', random_state=0, n_init=1, verbose=1)
prediction = kmeans.fit_predict(coverages_df.loc[:, ['coverage','GC_percentage']])


plt.scatter(x=coverages_df['GC_percentage'], y=coverages_df['coverage'], c=prediction)

plt.show()
# gc_content = pd.read_csv('longest_read_every_cluster_gc_content.tsv', delimiter='\t',header=None)
# df = gc_content[gc_content[2] < 30]
# print(df.describe())
#gc_content.hist(bins=2)
#plt.show()

# inicio = time.time()
# file = pd.read_csv("overlaps_talaro_18_07_2023.paf",delimiter='\t',header=None)
# file [10] = file[10].astype(int)
# file = file[file[10]>500]
# file_grouped = (file.groupby(by=0))#.size()
# print(file_grouped.size().sort_values(ascending=False)[:20].describe())
# print(file_grouped.ngroups)
# final = time.time()
# tiempo = final - inicio
# print("pandas", tiempo)
