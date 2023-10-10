from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import numpy as np
import pandas as pd


def set_minimun_cov(clusters_info: pd.DataFrame, coverage: str) -> int:
    ## Get minimum coverage
    kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(
        clusters_info["coverage"].array.reshape(-1, 1)
    )
    cov_fdp = kde.score_samples(clusters_info["coverage"].array.reshape(-1, 1))
    local_min = argrelextrema(cov_fdp, np.less)[0]

    ## Check if the user set a coverage
    if int(coverage) == -1:
        if len(local_min) < 1:
            min_coverage = 2
        else:
            min_coverage = clusters_info.iloc[max(local_min), :]["coverage"]
    else:
        min_coverage = coverage
    print("Covertura minima admitida: ", min_coverage)
    return min_coverage
