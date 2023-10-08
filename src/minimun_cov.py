from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema

def set_minimun_cov (clusters_info, args):
    ## Get minimum coverage
    kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(
        clusters_info["coverage"].array.reshape(-1, 1)
    )
    cov_fdp = kde.score_samples(clusters_info["coverage"].array.reshape(-1, 1))
    local_min = argrelextrema(cov_fdp, np.less)[0]
    
    ## Check if the user set a coverage
    if int(args[3]) > -1:
        if len(local_min) < 1:
            min_coverage = 2
        else:
            min_coverage = clusters_info.iloc[max(local_min), :]["coverage"]
    else:
        min_coverage = args[3]
    print("Covertura minima admitida: ", min_coverage)
    return min_coverage