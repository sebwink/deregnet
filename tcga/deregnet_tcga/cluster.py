from sklearn.cluster import k_means
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class KMeansClustering:
    def __init__(self, X, n_clusters):
        self.X = X.as_matrix()
        self.centroid, self.label, self.inertia = k_means(X.as_matrix(), n_clusters)
        self.samples = list(X.columns)
        self.clusters = [[gi for i, gi in enumerate(self.samples) if self.label[i] == k] for k in range(n_clusters)]
        self.cluster_indices = [[i for i, _ in enumerate(self.samples) if self.label[i] == k] for k in range(n_clusters)]
        self.index2sample = {i: sample for i, sample in enumerate(self.samples)}
        
    def __getitem__(self, i):
        if (not isinstance(i, int)) or i >= len(self.clusters):
            raise IndexError
        return self.clusters[i]
    
    def __iter__(self):
        for cluster in self.cluster:
            yield cluster
    
    def get_indices(self, i):
        if not isinstance(i, int) or i >= len(self.clusters):
            raise IndexError
        return self.cluster_indices[i]
    
    def plot(self, cluster_order=None, figsize=(12,10), **kwargs):
        # strange behaviour: only use once!
        def flatten_cluster_indices(inds):
            cli = inds[0]
            for clik in inds[1:]:
                cli += clik
            return cli
        cluster_order = cluster_order if cluster_order is not None else list(range(len(self.clusters)))
        cluster_indices = flatten_cluster_indices([self.get_indices(i) for i in cluster_order])
        reordering = np.array([len(cluster_indices) * [i] for i in cluster_indices])
        samples = [self.index2sample[i] for i in cluster_indices]
        _, ax = plt.subplots(figsize=figsize)
        sns.heatmap(pd.DataFrame(data=self.X[reordering, reordering.transpose()],
                                 index=samples,
                                 columns=samples),
                    ax=ax,
                    **kwargs)
        
    def cluster_sizes(self):
        return [len(cluster) for cluster in self.clusters]
