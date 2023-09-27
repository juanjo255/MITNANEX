from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

## I need to separate clusters of reads by coverage
## NOTE: so far only no plants organisms will be tackle.

def cluster_kmer_profiles (kmer_reduction:pd.DataFrame) -> np.ndarray :
    
    ## Clustering reads using K-means expecting 2 clusters ##
    kmeans = KMeans(
        n_clusters=2,
        max_iter=100,
        init="k-means++",
        random_state=0,
        n_init=1,
        verbose=1,
    )
    mt_prediction = kmeans.fit_predict(kmer_reduction[["comp1", "comp2"]])
    return mt_prediction