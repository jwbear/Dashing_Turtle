import sys, re
import platform
import os.path
import traceback
import numpy as np
import pandas as pd
import math
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import DashML.data_fx as dfx
import varnaapi


#### Must Load Varna Jar locally #####
varnaapi.set_VARNA(sys.path[1] + "/DashML/Deconvolution/Varna/VARNAv3-93.jar")

#TODO add lid to cluster images
#TODO convert to db
#TODO update save path file

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

save_path = None


# includes low/high confidence shapemap
# TODO: auto varna plots
def get_shapemap():
    df = dfx.get_shapemap()
    #inlier
    df['Predicted_Shape_Map'].loc[df['Reactivity_profile'] < .4] = 0
    df['Predicted_Shape_Map'].loc[df['Reactivity_profile'] >= .4] = 1
    df['Predicted_Shape_Map'].loc[df['Reactivity_profile'] >= .8] = 2
    df['Position'] = df['Position'] + 1
    sequences = df['Sequence_Name'].unique()
    df.to_csv('/Users/timshel/NanoporeAnalysis/DashML/ShapeMap/varna_maps.csv')
    for seq in sequences:
        dft = df[df['Sequence_Name'] == seq]
        dft = dft[['Position', 'Predicted_Shape_Map']]
        dft.to_csv('/Users/timshel/NanoporeAnalysis/DashML/VARNA/shapemap_' + seq + '.csv', index=False)

def get_avepredict():
    dfs = dfx.get_structure_ext()
    dfs.rename(columns={'Position': 'position', 'Sequence_Name': 'contig'}, inplace=True)
    dfs = dfs[['position', 'contig']]
    dfs['position'] = dfs['position'] + 1
    df = pd.read_csv('/Users/timshel/NanoporeAnalysis/DashML/Predictors/reactivity_ranking.csv')
    #inlier
    df['VARNA'].loc[df['Reactivity'] < 4] = 0
    df['VARNA'].loc[df['Reactivity'] >= 4] = 1
    df['VARNA'].loc[df['Reactivity'] >= 6] = 2
    df['position'] = df['position'] + 1
    df = dfs.merge(df, on=['position', 'contig'], how='left')
    df.fillna(0, inplace=True)
    sequences = df['contig'].unique()
    for seq in sequences:
        dft = df[df['contig'] == seq]
        dft = dft[['position', 'VARNA']]
        dft.to_csv('/Users/timshel/NanoporeAnalysis/DashML/VARNA/predict_' + seq + '.csv', index=False)

def plot_shapemap():
    #df = pd.read_csv("/Users/timshel/NanoporeAnalysis/DashML/ShapeMap/varna_RNAse_P.csv")
    #df['Position'] = df['Position'] + 1
    #varnaapi.load_config('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    # sequence = 'GUUAAUCAUGCUCGGGUAAUCGCUGCGGCCGGUUUCGGCCGUAGAGGAAAGUCCAUGCUCGCACGGUGCUGAGAUGCCCGUAGUGUUCGUGCCUAGCGAAUCCAUAAGCUAGGGCAGCCUGGCUUCGGCUGGGCUGACGGCGGGGAAAGAACCUACGUCCGGCUGGGAUAUGGUUCGAUUACCCUGAAAGUGCCACAGUGACGGAGCUCUAAGGGAAACCUUAGAGGUGGAACGCGGUAAACCCCACGAGCGAGAAACCCAAAUGAUGGUAGGGGCACCUUCCCGAAGGAAAUGAACGGAGGGAAGGACAGGCGGCGCAUGCAGCCUGUAGAUAGAUGAUUACCGCCGGAGUACGAGGCGCAAAGCCGCUUGCAGUACGAAGGUACAGAACAUGGCUUAUAGAGCAUGAUUAACGUC'
    # ss = "(((((((((((((((((((((.(((((((((....)))))))))...[.[[.[[[[[(((((((((.(((...).)))))).....((((((((((((........)))))))((((((((((....))))))))).)((.(((((((((((((..((((.....))))).))))))).))))).....(((((............((((((((....)))))))).........)))..))))))))))))))...((((.........))))...((((((((((((.(...)))))))))))))(((((((........)))))))........)))))))((((((((((..(((.....))).......))))))))))......]]]]]]]].).)))))))))))))..."
    # v = varnaapi.Structure(structure=ss, sequence=sequence)
    v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    v.savefig("RNAse_Pp.png")

 #get secondary structure for cluster
 # TODO try save yml
def save_bpseq(seq, cluster, structure):
    sequence = structure.sequence
    structure = structure.structure
    seqlen = lib.get_seqlen(seq)

    # create sequence save directory
    save_path_dir = save_path + seq
    if not os.path.exists(save_path_dir):
        os.makedirs(save_path_dir)
    save_path_dir = save_path_dir + '/'

    base_pairs = []
    left = []
    for i in range(len(structure)):
        if re.search('\(', structure[i]):
            left.append(i + 1)
        if re.search('\)', structure[i]):
            base_pairs.append([left.pop(), i + 1])
    base_pairs.sort(key=lambda tup: (tup[0], tup[1]))
    print(base_pairs)

    f = open(save_path_dir + seq + "_" + str(cluster) + ".bpseq", "w")
    f.write('# ' + seq + ' ' + str(cluster) + '\n')
    n = 0
    for i in range(1,seqlen+1):
        #print(i)
        if (n < len(base_pairs)) and (base_pairs[n][0] == i):
            f.write(str(i) + ' ' + sequence[i-1] + ' ' + str(base_pairs[n][1]) + '\n')
            n = n + 1
        else:
            bp2 = [x[1] for x in base_pairs]
            try:
                m = bp2.index(i)
                f.write(str(i) + ' ' + sequence[i-1] + ' ' + str(base_pairs[m][0]) + '\n')
            except ValueError:
                f.write(str(i) + ' ' + sequence[i-1] + ' ' + str(0) + '\n')
    f.write('\n')
    f.close()
    return
def get_vplot(seq):
    # get sequences and lengths
    seq_name = seq
    seq_len = lib.get_seqlen(seq)
    sequence = lib.get_sequence(seq)

    # create sequence save directory
    save_path_dir = save_path + seq
    if not os.path.exists(save_path_dir):
        os.makedirs(save_path_dir)
    save_path_dir = save_path_dir + '/'

    #get secondary structure for cluster
    def get_ss(cluster):
        ss = ""
        f = open(putative_structures +"/" + seq +"/" + seq + "_" + str(cluster) + ".out", "r")
        lines = f.readlines()
        for line in lines:
            if "MEA" in line:
                l = line.split(' ')
                ss = l[0]
                break
        f.close()
        return ss

    # dataframes of centroids
    df = pd.read_csv(data_path + seq + "_centroids.csv")
    position = np.array(np.arange(1, seq_len + 1, dtype=int))
    df['position'] = 0
    clusters = df['cluster'].unique()
    for cluster in clusters:
        print(cluster)
        df.loc[df['cluster'] == cluster, 'position'] = position
        ss = get_ss(cluster)
        v = varnaapi.Structure(structure=ss, sequence=sequence)
        v.update(resolution=10, zoom=1)
        save_bpseq(seq, cluster, v)
        out_fig = save_path_dir + seq + "-" + str(cluster) + ".png"
        #annotating high reactivity regions.
        #v.add_highlight_region(11, 21)
        #v.add_colormap(values=np.arange(1, 10), vmin=30, vmax=40, style='bw')
        #values is an array where each position indicates color 0-n
        #overall style is applied
        # annotating interactions
        #v.add_colormap(values=[2,5,5,5,5, 0, 0, 0, 0, 3,3 ,3],style='energy')
        #v.add_aux_BP(1, 10, color='red')
        v.savefig(out_fig)
        v.show()
        #sys.exit(0)


    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #df = pd.read_csv("/Users/timshel/NanoporeAnalysis/DashML/ShapeMap/varna_RNAse_P.csv")
    #df['Position'] = df['Position'] + 1
    #varnaapi.load_config('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    # sequence = 'GUUAAUCAUGCUCGGGUAAUCGCUGCGGCCGGUUUCGGCCGUAGAGGAAAGUCCAUGCUCGCACGGUGCUGAGAUGCCCGUAGUGUUCGUGCCUAGCGAAUCCAUAAGCUAGGGCAGCCUGGCUUCGGCUGGGCUGACGGCGGGGAAAGAACCUACGUCCGGCUGGGAUAUGGUUCGAUUACCCUGAAAGUGCCACAGUGACGGAGCUCUAAGGGAAACCUUAGAGGUGGAACGCGGUAAACCCCACGAGCGAGAAACCCAAAUGAUGGUAGGGGCACCUUCCCGAAGGAAAUGAACGGAGGGAAGGACAGGCGGCGCAUGCAGCCUGUAGAUAGAUGAUUACCGCCGGAGUACGAGGCGCAAAGCCGCUUGCAGUACGAAGGUACAGAACAUGGCUUAUAGAGCAUGAUUAACGUC'
    # ss = "(((((((((((((((((((((.(((((((((....)))))))))...[.[[.[[[[[(((((((((.(((...).)))))).....((((((((((((........)))))))((((((((((....))))))))).)((.(((((((((((((..((((.....))))).))))))).))))).....(((((............((((((((....)))))))).........)))..))))))))))))))...((((.........))))...((((((((((((.(...)))))))))))))(((((((........)))))))........)))))))((((((((((..(((.....))).......))))))))))......]]]]]]]].).)))))))))))))..."
    # v = varnaapi.Structure(structure=ss, sequence=sequence)
    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #v.savefig("RNAse_Pp.png")




#get_avepredict()
#get_shapemap()
#plot_shapemap()
# plot putative structures
get_vplot('HCV')
sys.exit(0)
