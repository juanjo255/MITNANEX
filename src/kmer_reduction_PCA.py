from sklearn.decomposition import PCA
import pandas as pd


def kmer_reduction(kmer_profiles_df:pd.DataFrame, clusters_info:pd.DataFrame, n_comp:int) -> pd.DataFrame:
    """Reduce dimensionaly using PCA

    Args:
        kmer_profiles_df (pd.DataFrame): Contains the kmer composition from the clusters
        clusters_info (pd.DataFrame): Contains the clusters
        n_comp (int): number of components

    Returns:
        pd.DataFrame: cluster_info columns plus the new dimensions
    """
    pca = PCA(n_components=n_comp)
    kmer_reduction_df = pd.DataFrame(
        pca.fit_transform(kmer_profiles_df.iloc[:, :-1]), columns=["comp1", "comp2"]
    )
    
    ## Merging the dataframe with ids and other relevant information ##
    kmer_reduction_df["ids"] = kmer_profiles_df["ids"]
    kmer_reduction_df = kmer_reduction_df.merge(
        clusters_info, how="left", left_on="ids", right_on="id_longest_read"
    )
    kmer_reduction_df.drop(columns="id_longest_read", inplace=True)
    return kmer_reduction_df
