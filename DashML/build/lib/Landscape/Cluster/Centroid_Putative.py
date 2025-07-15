import sys, re
import platform
import os.path
import traceback
import numpy as np
import pandas as pd
import varnaapi
from importlib.resources import files
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

######### Draws Varna Images of Predicted Secondary Structures #####
# TODO: Different from Native Images

#### Must Load Varna Jar locally #####
varna_path = files("DashML.Varna") / "VARNAv3-93.jar"
varnaapi.set_VARNA(varna_path)

#TODO add lid to cluster images
#TODO update save path file

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

save_path = '.'


 #get secondary structure for centroids
def save_bpseq(seq, seqlen, cluster, structure, out):
    sequence = structure.sequence
    structure = structure.structure

    base_pairs = []
    left = []
    for i in range(len(structure)):
        if re.search('\(', structure[i]):
            left.append(i + 1)
        if re.search('\)', structure[i]):
            base_pairs.append([left.pop(), i + 1])
    base_pairs.sort(key=lambda tup: (tup[0], tup[1]))
    #print(base_pairs)

    f = open(out + seq + "_" + str(cluster) + ".bpseq", "w")
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

##### Plot Predicted Secondary Structure of Centroids ####
def get_vplot(lids):
    df = dbsel.select_putativestructures(lids)
    # get sequences and lengths
    seq = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['contig'].unique()))
    seqlen = int(df['sequence_len'].unique()[0])
    sequence = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['sequence'].unique()))

    def draw_structure(df, method='hamming'):
        out = save_path + seq +'_'+ str(lids) + '/' + method + '/'
        if not os.path.exists(out):
            os.makedirs(out)
        print(out)
        clusters = df['cluster'].unique()
        for cluster in clusters:
            print(cluster)
            ss = str(df.loc[df['cluster'] == cluster, 'secondary'].unique()[0])
            v = varnaapi.Structure(structure=ss, sequence=sequence)
            v.update(resolution=10, zoom=1)
            out_fig = out + seq + "_" + str(cluster) + ".png"
            save_bpseq(seq, seqlen, cluster, v, out)

            # annotating high reactivity regions.
            # v.add_highlight_region(11, 21)
            # v.add_colormap(values=np.arange(1, 10), vmin=30, vmax=40, style='bw')
            # values is an array where each position indicates color 0-n
            # overall style is applied
            # annotating interactions
            # v.add_colormap(values=[2,5,5,5,5, 0, 0, 0, 0, 3,3 ,3],style='energy')
            # v.add_aux_BP(1, 10, color='red')
            v.savefig(out_fig)
            v.show()
            # sys.exit(0)

    # dataframes of centroids
    df_hamming = df.loc[df['method'] == 'hamming']
    if len(df_hamming) > 0:
        draw_structure(df_hamming, method='hamming')
    df_kmeans = df.loc[df['method'] == 'kmeans']
    if len(df_kmeans) > 0:
        draw_structure(df_kmeans, method='kmeans')

    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #df = pd.read_csv("/Users/timshel/NanoporeAnalysis/DashML/ShapeMap/varna_RNAse_P.csv")
    #df['Position'] = df['Position'] + 1
    #varnaapi.load_config('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    # sequence = 'GUUAAUCAUGCUCGGGUAAUCGCUGCGGCCGGUUUCGGCCGUAGAGGAAAGUCCAUGCUCGCACGGUGCUGAGAUGCCCGUAGUGUUCGUGCCUAGCGAAUCCAUAAGCUAGGGCAGCCUGGCUUCGGCUGGGCUGACGGCGGGGAAAGAACCUACGUCCGGCUGGGAUAUGGUUCGAUUACCCUGAAAGUGCCACAGUGACGGAGCUCUAAGGGAAACCUUAGAGGUGGAACGCGGUAAACCCCACGAGCGAGAAACCCAAAUGAUGGUAGGGGCACCUUCCCGAAGGAAAUGAACGGAGGGAAGGACAGGCGGCGCAUGCAGCCUGUAGAUAGAUGAUUACCGCCGGAGUACGAGGCGCAAAGCCGCUUGCAGUACGAAGGUACAGAACAUGGCUUAUAGAGCAUGAUUAACGUC'
    # ss = "(((((((((((((((((((((.(((((((((....)))))))))...[.[[.[[[[[(((((((((.(((...).)))))).....((((((((((((........)))))))((((((((((....))))))))).)((.(((((((((((((..((((.....))))).))))))).))))).....(((((............((((((((....)))))))).........)))..))))))))))))))...((((.........))))...((((((((((((.(...)))))))))))))(((((((........)))))))........)))))))((((((((((..(((.....))).......))))))))))......]]]]]]]].).)))))))))))))..."
    # v = varnaapi.Structure(structure=ss, sequence=sequence)
    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #v.savefig("RNAse_Pp.png")





# plot putative structures
#get_vplot('37')
#sys.exit(0)
