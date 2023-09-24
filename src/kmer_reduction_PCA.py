from sklearn.decomposition import PCA
import pandas as pd


def kmer_reduction(kmer_profiles_df:pd.DataFrame, clusters_info:pd.DataFrame) -> pd.DataFrame:
    pca = PCA(n_components=2)
    kmer_reduction_df = pd.DataFrame(
        pca.fit_transform(kmer_profiles_df.iloc[:, :-2]), columns=["comp1", "comp2"]
    )
    
    ## Merging the dataframe with ids and other relevant information ##
    kmer_reduction_df["ids"] = kmer_profiles_df["ids"]
    kmer_reduction_df = kmer_reduction.merge(
        clusters_info, how="left", left_on="ids", right_on="id_longest_read"
    )
    kmer_reduction.drop(columns="id_longest_read", inplace=True)
    return kmer_reduction_df
