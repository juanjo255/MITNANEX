from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

## I need to separate clusters of reads by coverage
## NOTE: so far only no plants organisms will be tackle.

def cluster_kmer_profiles (kmer_reduction_df:pd.DataFrame, n_clusters:int, max_iter:int) -> np.ndarray :
    
    ## Clustering reads using K-means expecting 2 clusters ##
    kmeans = KMeans(
        n_clusters=n_clusters,
        max_iter=max_iter,
        init="k-means++",
        random_state=0,
        n_init=1,
        verbose=0,
    )
    mt_prediction = kmeans.fit_predict(kmer_reduction_df[["comp1", "comp2"]])
    return mt_prediction