import sys, re
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering, KMeans
from scipy.spatial.distance import squareform
from scipy.spatial import distance
from kmodes.kmodes import KModes


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if platform.system() == 'Linux':
    ##### server #####
    data_path = "/home/jwbear/projects/def-jeromew/jwbear/dendrogram/Out/"
    save_path = "/home/jwbear/projects/def-jeromew/jwbear/dendrogram/Dendrogram/Dendrogram_Out/"
else:
    data_path = sys.path[1] + "/DashML/Deconvolution/Dendrogram/Putative_Structures/"
    out_path = sys.path[1] + "/DashML/Deconvolution/Out/"
    save_path = sys.path[1] + "/DashML/Deconvolution/Dendrogram/Putative_Structures/"


def get_mfe_distances(seq='HCV'):
    # native conformation mfes from RNAeval
    df_native = pd.read_csv(data_path + "native_mfes.txt", names=['sequence', 'mfe'],
                              dtype={'sequence':str, 'mfe': float},
                              skipinitialspace=True)
    # putative conformations mfes from RNAeval, multiple can take average?
    df_putative = pd.read_csv(data_path + "putative_mfes.txt", names=['sequence', 'cluster_type', 'cluster_number', 'mfe'],
                              dtype={'sequence':str, 'cluster_type': str, 'cluster_number':int, 'mfe': float},
                              skipinitialspace=True)
    #cluster sizes
    df_counts = pd.read_csv("/Users/timshel/structure_landscapes/DashML/Deconvolution/Dendrogram/Clusters/"
                            "cluster_counts.csv", names=['sequence', 'cluster_type', 'cluster_number', 'cluster_size'],
                            dtype={'sequence': str, 'cluster_type': str, 'cluster_number': int, 'cluster_size': int},
                            skipinitialspace=True)

    #merge data
    df_counts['mfe'] = 0
    df_putative['cluster_size'] = 0
    df = df_putative.merge(df_counts, on=['sequence', 'cluster_type', 'cluster_number'], how='left')
    df.drop(columns=['mfe_y', 'cluster_size_x'], inplace=True)
    df.rename(columns={'mfe_x':'mfe', 'cluster_size_y': 'cluster_size'}, inplace=True)
    df = df.sort_values(by=['sequence', 'cluster_type', 'cluster_number'])

    # get delta(native mfe, cluster mfe)
    # todo add seq variable loop
    native_mfe = df_native.loc[df_native['sequence']==seq, 'mfe'][0]
    df['native_mfe'] = native_mfe
    # negative means mfe of cluster is more stable than native conformation
    # a more negative mfe is more stable
    # a more positive mfe is less stable
    # native should be most negative, more negative the better in clusters
    df['delta_native'] = np.subtract(np.abs(df['native_mfe']),np.abs(df['mfe']))
    print(df)

#get_mfe_distances()
#sys.exit(0)



#get_cluster_hamming_distance()
#distance_matrix()
