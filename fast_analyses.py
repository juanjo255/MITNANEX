import time
from src.mitnanex import run
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from Bio import SeqIO
from sourmash import MinHash

## FASTA file is lighter
# convert_fq_to_fa(
#     "test/aedes_vexans_all_reads_subsample_sorted_length.fastq",
#     "test/aedes_vexans_all_reads_subsample_sorted_length.fasta",
# )
records_dict = SeqIO.to_dict(
        SeqIO.parse("test/aedes_vexans_all_reads_subsample_sorted_length.fasta", "fasta")
    )


inicio = time.time()
a = run()

# I need to plot the coverage of clusters
coverages_df = pd.DataFrame(
    {
        "coverage": [i.coverage for i in a.clusters],
        "longest_read_len": [i.longest_read_length for i in a.clusters],
        "id_longest_read": [i.longest_read_id for i in a.clusters],
        "id_cluster":[i.id_cluster for i in a.clusters],
    }
)

# Test on the first 30 clusters
prueba = coverages_df.sort_values(by='coverage', ascending=False).head(30)
for id in prueba['id_longest_read']:
    records_dict[id].seq

final = time.time()
tiempo = final - inicio
print("my code", tiempo)


## First view before Kmeans
# plt.scatter(x=coverages_df['GC_percentage'], y=coverages_df['coverage'])

# kmeans = KMeans(
#     n_clusters=3, max_iter=100, init="k-means++", random_state=0, n_init=1, verbose=1
# )
# prediction = kmeans.fit_predict(coverages_df.loc[:, ["coverage", "GC_percentage"]])


# plt.scatter(x=coverages_df["GC_percentage"], y=coverages_df["coverage"], c=prediction)

# plt.show()

# USING PANDAS
# file = pd.read_csv("test/overlaps_talaro_18_07_2023_sorted_containment.paf",delimiter='\t',header=None)
# file [10] = file[10].astype(int)
# print(file.describe())
# file = file[file[10]>500]
# print(file.describe())
# file_grouped = (file.groupby(by=0))
# print(file_grouped.head())
# print(file.sort_values(by=1, ascending=False).head())
