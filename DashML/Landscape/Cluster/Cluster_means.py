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
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import distance
import DashML.Database_fx.Insert_DB as dbins


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
def distancematrix(df, seq='hcv'):
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

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

#### method ultimately decided by which method most closely approximates the desired function
def get_clusters(X, seq, read_list, tit, clust_num):
    # Heirarchical Clustering
    #perform clustering
    model = AgglomerativeClustering(n_clusters=clust_num, linkage='complete', metric='precomputed',
                                    compute_distances=True, compute_full_tree=True)
    model = model.fit(X)
    df = pd.DataFrame(columns=['read_id', 'cluster'])
    df['read_id'] = read_list
    df['cluster'] = model.labels_
    # plot the top three levels of the dendrogram
    plt.figure(figsize=(10, 7))
    plt.title("KMeans Hierarchical Clustering Dendrogram")
    plt.xlabel("Number of points in node (or index of point).")
    ax = plt.gca()
    plot_dendrogram(model, truncate_mode="level", p=5, ax=ax)
    img_path = save_path + seq + tit + "_dendrogram.png"
    plt.savefig(img_path, dpi=600)
    #plt.showblock=False)
    return df, img_path


# hybrid correlation matrix using distance instance of correlation
def clustermap(X, seq='HCV', tit='Kmeans'):
    # cannot contain nan
    X = np.nan_to_num(X)
    z_col = linkage(np.transpose(squareform(X)), method='complete')
    print(z_col.shape)
    z_row = linkage(squareform(X), method='complete')
    print(z_row.shape)
    g = sns.clustermap(data=X, row_linkage=z_row, col_linkage=z_col, figsize = (8,8))
    #g.ax_heatmap.set_xlabel("Cluster Labels", labelpad=1, fontweight='extra bold')
    #g.ax_heatmap.set_ylabel("Clusters Labels", labelpad=1, fontweight='extra bold')
    #plt.title(tit + " Inter-Cluster Distances ")
    plt.ylabel("Distance (Min-Max)")
    g.savefig(save_path + seq + tit + '_clustermap.png', dpi=600)
    #plt.showblock=False)
    #print(dir(g))
    #print(g.data2d)
    return


def corrmatrix(seq='HCV', data=None, tit="Kmeans"):
    g = sns.heatmap(np.corrcoef(data))
    g.set_xlabel("Reads")
    g.set_ylabel("Reads")
    g.set_title( "Kmeans Read Correlations")
    g.get_figure().savefig(save_path + seq + tit + '_read_corr_matrix.png', dpi=600)
    img_path = save_path + seq + tit + '_read_corr_matrix.png'
    #plt.showblock=False)
    #print(np.corrcoef(data))
    return img_path

#### kmeans to get centroids ####
def get_centroids(X, clusters, seq):
    # set cluster labels to reads
    df = X.merge(clusters, on=['read_id'], how='left')
    df['method'] = 'kmeans'
    seq_ids = df['LID'].unique().flatten()
    # save all clusters
    dbins.insert_clusters(df, seq_ids)
    df_centroid = pd.DataFrame(columns=['LID', 'cluster', 'position', 'centroid'])
    seqlen = int(df['position'].max())
    for i in clusters['cluster'].unique():
        dt = df.loc[df['cluster']==i]
        lids = dt['LID'].unique()
        dt = dt[['read_index', 'position', 'Reactivity']]
        dt = dt.sort_values(by=['read_index', 'position'])
        dt = (dt.pivot(index=["read_index"], columns="position", values="Reactivity")
               .rename_axis(columns=None)
               .reset_index())
        dt.fillna(value=0, inplace=True)
        dt.drop(columns=['read_index'], inplace=True)
        X = dt.to_numpy()
        # get kmeans for all reads in cluster
        km = KMeans(n_clusters=1, random_state=0, n_init="auto").fit(X)
        for lid in lids:
            df_centroid = pd.concat([df_centroid, pd.DataFrame({'LID': lid, 'cluster': i, 'position': np.arange(0, seqlen+1),
                                                            'centroid': np.array(km.cluster_centers_).flatten()})])
    df_centroid['method'] = 'kmeans'
    df_centroid['contig'] = re.sub('\[|\]|\\s|\\n+|\'|\"', '', seq)

    # save centroids
    dbins.insert_centroids(df_centroid,seq_ids)
    return

###### Heatmap of Reads Against Summed Reactivity Across different Metrics ######
def heatmap(df, seq, tit):
    fig = plt.gcf()
    plt.figure(figsize=(60, 60))
    ax = plt.gca()
    df = df[['position', 'contig', 'read_id', 'Reactivity']]
    result = df.pivot(index="position", columns="read_id", values="Reactivity")
    seq_len = len(df['position'].unique())
    sns.heatmap(result, cmap='crest',
                      yticklabels=range(0,seq_len), ax=ax)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.xlabel("Read Number", fontsize=100, weight='bold')
    plt.ylabel("Position", fontsize=100, weight='bold')
    plt.title(seq + " Reactivity", fontsize=100, weight='bold')
    img_path = save_path + seq + tit + '_heatmap.png'
    plt.savefig(img_path, dpi=300)
    #plt.showblock=False)
    return img_path


def den_array(df, clust_num, plot=True):
    tit = '_kmeans_' + str(df['LID'].unique())
    seq = str(df['contig'].unique())
    seq = re.sub('\[|\]|\\s|\\n+|\'|\"', '', seq)
    distk, reads, read_list = distancematrix(df, seq=seq)
    clusters, img_path_dend = get_clusters(X=distk, read_list=read_list, seq=seq, tit=tit, clust_num=clust_num)
    #centroids
    get_centroids(X=df, clusters=clusters, seq=seq)
    if plot:
        img_path_heatmap = heatmap(df, seq, tit)
        #get ordered similarity from most to least similar for clusters
        clustermap(distk, seq=seq, tit=tit)
        # correlation between reads
        img_path_corr = corrmatrix(seq=seq, data=reads, tit=tit)
        return img_path_dend, img_path_heatmap, img_path_corr
    del df
    return None
