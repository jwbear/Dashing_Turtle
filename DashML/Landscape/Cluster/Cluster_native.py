import sys, re
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from sklearn.metrics.cluster import adjusted_mutual_info_score
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

matplotlib.use("Agg")
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

save_path = None

#reactivity format for RNAcofold
def scale_reactivities(reactivities):
    min = reactivities.min()
    max = reactivities.max()
    smin = 0
    smax = 2

    reactivities = ((reactivities - min) / (max - min)) * (smax - smin) + smin
    return reactivities


# Create distance matrices between native sequence and cluster centroids using hamming and kmeans

def distance_matrix(df_kmeans, df_hamming, native_lid, seq="HCV"):
    native_rx = get_native_rx(native_lid)

    ##### Compute Distance Between centroids and native #####
    #reactivity kmeans distance
    dx = df_kmeans.sort_values(by=['LID','cluster', 'position'])
    dx = (dx.pivot(index=["cluster"], columns="position", values="reactivity")
         .rename_axis(columns=None)
         .reset_index())
    dx.fillna(value=0, inplace=True)
    dx.drop(columns=['cluster'], inplace=True)
    mm = dx.to_numpy()

    # kmeans distance matrix
    km = distance.cdist(native_rx, mm, 'euclidean')
    k_corr = mm

    #save centroid distances
    df = df_kmeans.sort_values(by=['LID', 'cluster','position'])
    df = df[['LID', 'cluster']].reset_index(drop=True)
    df = df.groupby(by=['LID', 'cluster']).mean().reset_index()
    df['distance'] = np.array(km).flatten()
    df['method'] = 'kmeans'
    dbins.insert_centroid_distance(df)

    # reactivity kmode distance
    dx = df_hamming.sort_values(by=['LID', 'cluster', 'position'])
    dx['reactivity'] = np.where(dx['reactivity'] == -1, 1, 0)
    dx = (dx.pivot(index=["cluster"], columns="position", values="reactivity")
          .rename_axis(columns=None)
          .reset_index())
    dx.fillna(value=0, inplace=True)
    dx.drop(columns=['cluster'], inplace=True)
    mm = dx.to_numpy()

    # hamming distance matrix
    ham = distance.cdist(native_rx, mm, 'hamming')
    ham_corr = mm

    # save centroid distances
    df = df_hamming.sort_values(by=['LID', 'cluster', 'position'])
    df = df[['LID', 'cluster']].reset_index(drop=True)
    df = df.groupby(by=['LID', 'cluster']).mean().reset_index()
    df['distance'] = np.array(ham).flatten()
    df['method'] = 'hamming'
    dbins.insert_centroid_distance(df)

    #correlation between clusters matrix
    ham_map = distance.cdist(mm, mm, 'hamming')

    return k_corr, ham_corr, km, ham

# get reactivities for native structures using hamming binary reactivities
# remove sections with no data for each cluster, only compare non-zero positions
# TODO: get reactivities for native structures using native reactivities calculated in Native_Nanopore
# remove sections with no data for each cluster, only compare non-zero positions
def get_native_rx(lid):
    df_native = dbsel.select_structure(lid)
    df_native = df_native[['base_type']]
    df_native['base_type'] = np.where(df_native['base_type']=='S', 1, 0)
    mm = df_native.to_numpy()
    return mm.T

###### Correlation between Representatives Centroids ######
def corrmatrixh(seq='HCV', data=None, tit="Hamming", clust_num=5):
    plt.clf()
    plt.figure(figsize=(20, 20))
    print(data.shape)
    arr = np.zeros((clust_num, clust_num))
    for i in range(0, clust_num-1):
        for j in range(0, clust_num-1):
            arr[i,j] = adjusted_mutual_info_score(data[i,:], data[j,:])
    g = sns.heatmap(arr, annot=True)
    g.set_xlabel("Clusters")
    g.set_ylabel("Clusters")
    g.set_title( tit + " Cluster Correlations (Adjusted Mutual Information)")
    g.get_figure().savefig(save_path + seq +'_'+ tit + '_corr_matrix.png', dpi=600)
    #plt.showblock=False)
    #print(np.corrcoef(data))

def corrmatrixk(seq='HCV', data=None, tit="Kmeans", clust_num=5):
    plt.clf()
    plt.figure(figsize=(20, 20))
    arr = np.zeros((clust_num, clust_num))
    for i in range(0, clust_num - 1):
        for j in range(0, clust_num - 1):
            arr[i, j] = adjusted_mutual_info_score(data[i, :], data[j, :])
    #g = sns.heatmap(np.nan_to_num(np.corrcoef(data)), annot=True)
    g = sns.heatmap(arr, annot=True)
    g.set_xlabel("Clusters")
    g.set_ylabel("Clusters")
    g.set_title( tit + " Cluster Correlations (Adjusted Mutual Information)")
    img_path = save_path + seq +'_'+ tit + '_corr_matrix.png'
    g.get_figure().savefig(img_path, dpi=600)
    #plt.showblock=False)
    return img_path



###### Heatmap of Cluster Representatives Against native structure ######
def heatmap(X, seq="HCV", method="Kmeans"):
    fig = plt.gcf()
    plt.figure(figsize=(20, 8))
    plt.title( seq + " " + method + " Centroid Distance from Native Structure.\n")
    ax = sns.heatmap(X, cmap='crest', annot=True, square=True)
    ax.set(xlabel="Clusters", ylabel=seq)
    plt.savefig(save_path + seq + '_' + method + '_heatmap.png', dpi=600)
    #plt.show)


def den_array(lids, native_lid, plot=True):
    df = dbsel.select_centroids(lids)
    seq = str(df['contig'].unique())
    seq = re.sub('\[|\]|\\s|\\n+|\'|\"', '', seq)
    df_kmeans = df[df['method']=='kmeans']
    df_hamming = df[df['method']=='hamming']
    clust_num = len(df['cluster'].unique())

    k_corr, ham_corr, dist_k, dist_ham = distance_matrix(df_kmeans, df_hamming, native_lid=native_lid, seq=seq)
    if plot:
        # heatmap of cluster distances from native sequences
        heatmap(dist_ham, method="Hamming")
        heatmap(dist_k, method="Kmeans")
        # todo: get ordered similarity from most to least similar for clusters
        # TODO: correlation between reads or clusters???
        corrmatrixh(seq=seq, data=ham_corr, tit="Hamming", clust_num=clust_num)
        img_corr_means = corrmatrixk(seq=seq, data=k_corr, tit="Kmeans",  clust_num=clust_num)
        return img_corr_means
    return None
