import sys, re
import os.path
import platform
import traceback
import numpy as np
import pandas as pd
import math
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if platform.system() == 'Linux':
    ##### server #####
    data_path = "/home/jwbear/projects/def-jeromew/jwbear/Deconvolution/"
    save_path = "/home/jwbear/projects/def-jeromew/jwbear/Deconvolution/Out/"
else:
    data_path = sys.path[1] + "/DashML/Deconvolution/Dendrogram/"
    save_path = sys.path[1] + "/DashML/Deconvolution/Dendrogram/Tables/"
    save_fig = sys.path[1] + "/DashML/Deconvolution/Dendrogram/Figures/"


# creates a long table of putative mfe distances from native
def mfe_table(seq="HCV"):
    df_native = pd.read_csv(data_path + "Native/native_mfes.txt", names=['Sequence', 'MFE'])
    df_put = pd.read_csv(data_path + "Putative_Structures/" + seq + "_putative_mfes.txt", names=['Sequence', 'Method',
                                                                                     'Cluster', 'MFE'])
    for seq in df_native['Sequence'].unique():
        native_mfe = df_native.loc[df_native['Sequence']==seq, 'MFE'].to_numpy()
        print(native_mfe)
        for met in df_put['Method'].unique():
            dt = df_put.loc[(df_put['Sequence']==seq) & (df_put['Method']==met)]
            dt.drop(columns=['Sequence', 'Method'], inplace=True)
            dt.sort_values(by=['Cluster','MFE'])
            #todo euclidean distance
            dt['Distance'] = np.abs(dt['MFE']) - np.abs(native_mfe)
            latex = dt.to_latex(index=False, float_format="{:.2f}".format)
            latex = str(latex).replace("_", "\\_")
            #remove included header
            hindex = latex.find("\midrule")
            latex = latex[hindex:]
            header = ("\\begin{longtable}{|l|c|c|}\n" +
                      "\\caption*{\\textbf{" + seq +" "+ met + " Cluster Distance from Native MFE}}\n" +
                      "\\label{tab:" + seq + "nativemfedistances}\\\\\n")
            header = header + ("\\hline \n" +
                      "\\thead{Cluster} & \\thead{Putative MFE} & \\thead{Distance}\\\\"
                      "\n\hline\n")
            hindex = latex.find("\end{tabular}")
            latex = latex[:hindex]
            latex = header + latex
            latex = latex + "\\end{longtable}"
            with open(save_path + seq +" "+ met +  "_distances_from_native.tex", 'w+') as f:
                f.write(latex)
    return

#TODO: network x of clusters to show compact vs not compactness in mfes or structural distance????

# within paper display of minimum mfe distances and associated clusters
def mfe_table(seq="HCV"):
    df_native = pd.read_csv(data_path + "Native/native_mfes.txt", names=['Sequence', 'MFE'])
    df_put = pd.read_csv(data_path + "Putative_Structures/" + seq+ "/" + seq + "_putative_mfes.txt", names=['Sequence', 'Method',
                                                                                     'Cluster', 'MFE'])
    for seq in df_native['Sequence'].unique():
        native_mfe = df_native.loc[df_native['Sequence']==seq, 'MFE'].to_numpy()
        print(native_mfe)
        for met in df_put['Method'].unique():
            dt = df_put.loc[(df_put['Sequence']==seq) & (df_put['Method']==met)]
            dt.drop(columns=['Sequence', 'Method'], inplace=True)
            dt.sort_values(by=['Cluster','MFE'])
            dt['Distance'] = np.abs(dt['MFE']-native_mfe)
            latex = dt.to_latex(index=False, float_format="{:.2f}".format)
            latex = str(latex).replace("_", "\\_")
            #remove included header
            hindex = latex.find("\midrule")
            latex = latex[hindex:]
            header = ("\\begin{longtable}{|l|c|c|}\n" +
                      "\\caption*{\\textbf{" + seq +" "+ met + " Cluster Distance from Native MFE}}\n" +
                      "\\label{tab:" + seq + "nativemfedistances}\\\\\n")
            header = header + ("\\hline \n" +
                      "\\thead{Cluster} & \\thead{Putative MFE} & \\thead{Distance}\\\\"
                      "\n\hline\n")
            hindex = latex.find("\end{tabular}")
            latex = latex[:hindex]
            latex = header + latex
            latex = latex + "\\end{longtable}"
            with open(save_path + seq +" "+ met +  "_distances_from_native.tex", 'w+') as f:
                f.write(latex)
    return


def latex_table(df, seq, fname="_base_file", cols=None):
    #cols = ['Sequence', 'Canonical Modification Rates', 'Non-Canonical Modification Rates']
    #df = pd.DataFrame(columns=cols)
    latex = df.to_latex(index=False,float_format="{:.2f}".format)
    latex = str(latex).replace("_", "\\_")
    with open(save_path + seq + "_" + fname + ".tex", 'w+') as f:
        f.write(latex)


def mfe_cluster_graph(seq='HCV', method='kmeans'):
    df_native = pd.read_csv(data_path + "Native/native_mfes.txt", names=['sequence', 'MFE'])

    df_put = pd.read_csv(data_path + "Putative_Structures/" + seq + "/" + seq + "_putative_mfes.txt",
                         names=['sequence', 'method','cluster', 'MFE'])
    df_counts = pd.read_csv(data_path + "Clusters/" +seq + "_" +method+ "_cluster_counts.csv", names=['sequence', 'method',
                                                                                     'cluster', 'count'])
    df_counts = df_counts.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    df_counts['method'] = method
    df_counts['sequence'] = seq
    #distances from native structure
    df_struct = pd.read_csv(data_path + "Native/" +seq+ "_" +method+ "_structure_distance.csv", header=0)
    df_struct['method'] = method
    df_struct.rename(columns={'Distance': 'Structural Distance', 'Cluster':'cluster'}, inplace=True)
    df = df_put.merge(df_struct, on=['cluster', 'method'], how='inner')
    #distances from native mfe
    native_mfe = df_native.loc[df_native['sequence'] == seq, 'MFE'].to_numpy()
    df['MFE Distance'] = df['MFE'].apply(lambda x: int(np.abs(x - native_mfe)))
    df = df.merge(df_counts, on=['method', 'sequence', 'cluster'], how='left')
    sizes = df['count'].to_numpy()
    ax = df.plot.scatter(y='MFE Distance', x='Structural Distance',c='cluster', colormap = 'viridis_r', s=sizes)
    fig_name = seq + " " + method.title() + " Structure and MFE Distance"
    plt.title(fig_name)
    plt.savefig(save_fig+fig_name, dpi=600)
    #plt.show)

    ### make tables of min values ###
    ### average mfe, structural distance
    df = df.groupby(by=['sequence', 'method', 'cluster']).mean()
    df.sort_values(by=['Structural Distance', 'MFE Distance'], inplace=True)
    df.drop(columns=['MFE'], inplace=True)
    df.rename(columns={'count': 'Cluster Size'}, inplace=True)
    tot_clusters = df['Cluster Size'].sum()
    table_info = "Native MFE: " + str(native_mfe) + "\n" + \
        "Average MFE: " + str(df['MFE Distance'].mean()) + "\n" + \
        "Average Structural Distance: " + str(df['Structural Distance'].mean()) + "\n" + \
        "Total Clusters: " + str(tot_clusters)
    df['% Landscape'] = (df['Cluster Size']/tot_clusters) * 100
    print(df.sort_values(by=['% Landscape', 'Structural Distance', 'MFE Distance'], inplace=True, ascending=False))
    df = df.head(5)
    latex_table(df, seq, fname=method+"_"+"5", cols =df.columns)
    with open(save_path + seq + '_' + method + '_' + '5_table_info.txt', 'w') as f:
        f.write(table_info)
    print(df)

mfe_cluster_graph(seq='T_thermophila', method='kmeans')
mfe_cluster_graph(seq='T_thermophila', method='hamming')
mfe_table(seq='T_thermophila')
sys.exit(0)
