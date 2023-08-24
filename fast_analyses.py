import time
from src.mitnanex import run
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# inicio = time.time()
# a = run()

# # I need to plot the coverage of clusters 
# coverages_df = pd.DataFrame({"coverage":[i.coverage for i in a.clusters], "longest_read_len":[i.longest_read_length for i in a.clusters], "id_longest_read": [i.longest_read_id for i in a.clusters] })
# coverages_df.sort_values(by='coverage',inplace=True, ascending=False)
# print(coverages_df.duplicated('id_longest_read').sum())
# final = time.time()
# tiempo = final - inicio
#print("my code", tiempo)

gc_content = pd.read_csv('../longest_read_every_cluster_gc_content.tsv', delimiter='\t',header=None)
df = gc_content[gc_content[2] < 30]
print(df.describe())
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
