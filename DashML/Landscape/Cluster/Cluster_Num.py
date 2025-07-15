import os
import sys
import re
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import distance


#### Get Average Number of Optimal Clusters for Predictions ####


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

save_path=None

#reactivity format for RNAcofold
def scale_reactivities(reactivities):
    min = reactivities.min()
    max = reactivities.max()
    smin = 0
    smax = 2

    reactivities = ((reactivities - min) / (max - min)) * (smax - smin) + smin
    return reactivities

# Create linkage matrix and then plot the dendrogram
def distancematrix_modes(df, seq='hcv'):
    ##### Compute Distance Between Reads #####

    # reactivity kmode hamming
    dx = df.sort_values(by=['read_id', 'position'])
    read_list = dx['read_id'].unique()
    dx = (dx.pivot(index=["read_id"], columns="position", values="Predict")
          .rename_axis(columns=None)
          .reset_index())
    dx.fillna(value=0, inplace=True)
    dx = dx.drop(columns=['read_id'])
    mm = dx.to_numpy()

    #clustermap
    ham = distance.cdist(mm, mm, 'hamming')
    #ham = np.divide(ham,np.linalg.norm(ham))
    #print(ham)
    return ham, mm, read_list

def distancematrix_means(df, seq='hcv'):
    ##### Compute Distance Between Reads #####

    #reactivity kmeans distance
    dx = df.sort_values(by=['read_id', 'position'])
    read_list = dx['read_id'].unique()
    dx = (dx.pivot(index=["read_id"], columns="position", values="Reactivity")
         .rename_axis(columns=None)
         .reset_index())
    dx.fillna(value=0, inplace=True)
    dx = dx.drop(columns=['read_id'])

    ## create matrix for euclidean distance
    reads = dx.to_numpy()
    #mm = np.divide(mm, np.linalg.norm(mm))

    ## calculate distances
    distk = distance.cdist(reads, reads, 'euclidean')

    return distk, reads, read_list

def optimal_dendrogram(X, seq, tit, linkage_method='complete', threshold='auto', percent=.2):
    distance_matrix = pdist(X)
    Z = linkage(distance_matrix, method=linkage_method)

    plt.figure(figsize=(8, 5))
    dendrogram(Z)

    if 'kmeans' in tit:
        percent = .45

    if threshold == 'auto':
        max_height = np.max(Z[:, 2])
        threshold_val = percent * max_height
        num_clusters = np.sum(Z[:, 2] > threshold_val) + 1
        plt.axhline(y=threshold_val, color='r', linestyle='--',
                    label=f'{int(percent * 100)}% Height = {threshold_val:.2f}')
        plt.legend()
        print(f"Estimated number of clusters at {percent * 100}% height threshold: {num_clusters}")

    elif isinstance(threshold, (int, float)):
        num_clusters = np.sum(Z[:, 2] > threshold) + 1
        plt.axhline(y=threshold, color='r', linestyle='--', label=f'Threshold = {threshold}')
        plt.legend()
        print(f"Estimated number of clusters at threshold {threshold}: {num_clusters}")

    plt.title('Dendrogram')
    plt.xlabel('Read Index')
    plt.ylabel('Distance')
    plt.savefig(save_path + seq + tit + '_dendrogram.png', dpi=300)
    #plt.showblock=False)


def get_optimal_clusters(X, seq, read_list, tit, plot=False):
    # Heirarchical Clustering
    # perform clustering
    s_scores = []
    c_scores = []
    wss_values = []
    n_samples = len(X)
    max = min(50, n_samples)
    for k in np.arange(2,max, 1):
        model = AgglomerativeClustering(n_clusters=k, linkage='complete', metric='precomputed', compute_distances=True,
                                        compute_full_tree=True, distance_threshold=None)
        model = model.fit_predict(X)
        c_scores.append(1 + np.amax(model))
        #print(np.unique(model))
        #print(f"Number of clusters = {1 + np.amax(model)}")

        ### cluster evaluation
        # Silhouette Score maximal score
        if (1 + np.amax(model)) >= 2 and (k < len(np.unique(model)) < len(X)):
            silhouette = silhouette_score(X, model, metric='hamming')
            prev = silhouette
            s_scores.append(silhouette)
        else:
            s_scores.append(0)

         # Total Within-Cluster Sum of Squares (WSS)
        wss = np.sum([np.sum((X[model == c] - X[model == c].mean(axis=0)) ** 2) for c in np.unique(model)])
        wss_values.append(wss)
        #print('Silhouette Score:', silhouette, "WSS:", wss)

    print('Silhouette Score:', len(s_scores), "WSS:", len(wss_values))
    if plot:
        # Plot Elbow Curve
        plt.figure(figsize=(6, 4))
        plt.plot( c_scores, s_scores, marker='o', label='Silhouette Score')
        plt.title('Elbow Method using Silhouette Score')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Score')
        img_sil = save_path + seq + tit + '_Silhouette.png'
        plt.savefig(img_sil, dpi=300)
        #plt.showblock=False)

        # Plot Elbow Curve
        plt.figure(figsize=(6, 4))
        plt.plot( c_scores, wss_values, marker='*', label='WSS')
        plt.title('Elbow Method using WSS')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Score')
        img_wss = save_path + seq + tit + '_WSS.png'
        plt.savefig(img_wss, dpi=300)
        #plt.showblock=False)
    return img_sil, img_wss

def den_array(df, plot=True):
    seq = str(df['contig'].unique())
    seq = re.sub('\[|\]|\\s|\\n+|\'|\"', '', seq)
    ####### evaluate clusters

    ### modes
    tit = '_hamming_' + str(df['LID'].unique())
    distk, reads, read_list = distancematrix_modes(df, seq=seq)
    img_sil, img_wss = get_optimal_clusters(X=distk, read_list=read_list, seq=seq, tit=tit, plot=plot)
    optimal_dendrogram(X=distk, seq=seq, tit=tit)

    #### means
    tit = '_kmeans_' + str(df['LID'].unique())
    distk, reads, read_list = distancematrix_means(df, seq=seq)
    # evaluate clusters
    img_sil_k, img_wss_k = get_optimal_clusters(X=distk, read_list=read_list, seq=seq, tit=tit, plot=plot)
    optimal_dendrogram(X=distk, seq=seq, tit=tit)

    del df
    return img_sil_k, img_wss_k
