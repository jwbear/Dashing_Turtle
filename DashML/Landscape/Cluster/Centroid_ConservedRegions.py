#### get conserved regions for structural landscape #####
#### use centroids (kmeans/hamming) and count across landscape #####
#### set threshold percentage of landscape after calculation ####
#### conserved regions are more interesting for conformational or interaction changes ####

import os
import re
import platform
import sys
import pandas as pd
import numpy as np
import subprocess
from importlib.resources import files
import varnaapi
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

#### Must Load Varna Jar locally #####
varna_path = files("DashML.Varna") / "VARNAv3-93.jar"
varnaapi.set_VARNA(varna_path)


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
save_path='.'

def define_regions(df, cons_thresh,region_size, gap=1):
    cbp = []  # conserved bp region
    css = []  # conserved ss region

    #mask
    modsrx = ''.join(map(str,np.where(df['percent_mod']>=cons_thresh, 1, 0)))
    modsurx = ''.join(map(str, np.where(df['percent_mod'] < cons_thresh, 1, 0)))
    #set regions size
    regs = ""
    for i in range(0, region_size):
        regs = regs + "1"
    #evaluate regex
    mod_pattern = "("+ regs + "+)|(1.{1," + str(gap) + "}"+ (regs[0:len(regs)-2] + "+" if len(regs)>=3 else "")  +")"
    mod_pattern_urx = "(" + regs + "+)" #no gaps in urx, more likely to be a rx region
    #get regions, 0 indexed
    rx = re.finditer(mod_pattern, modsrx)
    urx = re.finditer(mod_pattern_urx, modsurx)
    for y in rx:
        i  = y.start()
        while i <= y.end():
            css.append(i)
            i= i + 1

    for y in urx:
        i  = y.start()
        while i <= y.end():
            cbp.append(i)
            i= i + 1
    return css, cbp

def insert_regions(df, css, cbp):
    df = df.groupby(['LID', 'contig', 'position', 'method']).size().reset_index()
    df= df[['LID', 'contig', 'position', 'method']]
    df['is_conserved'] = np.where(np.isin(df['position'].to_numpy(), cbp), 1, 0)
    df['is_conserved'] = np.where(np.isin(df['position'].to_numpy(), css), -1, df['is_conserved'])
    #insert
    dbins.insert_centroid_regions(df)
    return

def get_conserved_hamming(df, cons_thresh, region_size, gap):
    df_hamming = df.loc[df['method'] == 'hamming']
    clust_num = len(df_hamming['cluster'].unique())
    dt = (df_hamming.groupby(['position','method', 'centroid'], observed=False)['centroid'].size().
                       unstack(fill_value=0).reset_index())
    if -1 in dt.columns:
        dt['percent_mod'] = dt[-1]/clust_num
    else:
        dt['percent_mod'] = 0
    if 1 in dt.columns:
        dt['percent_unmod'] = dt[1]/clust_num
    else:
        dt['percent_unmod'] = 0
    css, cbp = define_regions(dt, cons_thresh=cons_thresh, region_size=region_size, gap=gap)
    insert_regions(df_hamming, css, cbp)
    return css, cbp

def get_conserved_kmeans(df, cons_thresh, region_size, gap, km_threshold=.4):
    #convert nonbinary df (kmeans)
    df.loc[(df['centroid'] < km_threshold) & (df['method'] == 'kmeans'), 'centroid'] = 1
    df.loc[(df['centroid'] >= km_threshold) & (df['method']=='kmeans') & (df['centroid'] !=1), 'centroid'] = -1
    df_km = df.loc[df['method'] == 'kmeans']
    clust_num = len(df_km['cluster'].unique())
    dt = (df_km.groupby(['position', 'method', 'centroid'], observed=False)['centroid'].size().
          unstack(fill_value=0).reset_index())
    if -1 in dt.columns:
        dt['percent_mod'] = dt[-1] / clust_num
    else:
        dt['percent_mod'] = 0
    if 1 in dt.columns:
        dt['percent_unmod'] = dt[1] / clust_num
    else:
        dt['percent_unmod'] = 0
    css, cbp = define_regions(dt, cons_thresh=cons_thresh, region_size=region_size, gap=gap)
    insert_regions(df_km, css, cbp)
    return css, cbp

#get secondary structure
def save_bpseq(lid, seq, seqlen, structure, metric):
    sequence = structure.sequence
    structure = structure.structure
    seqlen = seqlen
    f = None

    base_pairs = []
    left = []
    for i in range(len(structure)):
        if re.search('\(', structure[i]):
            left.append(i + 1)
        if re.search('\)', structure[i]):
            base_pairs.append([left.pop(), i + 1])
    base_pairs.sort(key=lambda tup: (tup[0], tup[1]))
    #print(base_pairs)
    try:
        f = open(save_path + str(lid) + "_" + seq + "_" + metric + ".bpseq", "w")
        f.write('# ' + seq + ' ' + metric + '\n')
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
    except Exception as e:
        raise Exception(str(e))
        sys.stderr.write(str(e))
        return None
    finally:
        if f is not None:
            f.close()
    return

def get_vplot(lid, seq, seq_len, sequence, ss, cons_ss, cons_bp, metric):
    v = varnaapi.Structure(structure=ss, sequence=sequence)
    v.update(resolution=10, zoom=1, algorithm='radiate', flat=True)
    save_bpseq(lid, seq, seqlen=seq_len, structure=v, metric=metric)
    out_fig = save_path + str(lid) + "_" + seq + "_" + metric + "_conserved.png"
    v.dump_param(save_path + str(lid) + "_" + seq + "_" + metric + "_conserved.yml")
    #annotating high reactivity regions.
    for r in cons_ss:
        if r < seq_len:
            r = r + 1
        v.add_highlight_region(r, r, fill='#f16849', outline='#f16849')
    #conserved inaccessible regions
    for r in cons_bp:
        if r < seq_len:
            r = r + 1
        v.add_highlight_region(r, r, fill='#c5def2', outline='#c5def2')
    #v.add_colormap(values=np.arange(1, 10), vmin=30, vmax=40, style='bw')
    #values is an array where each position indicates color 0-n
    #overall style is applied
    # annotating interactions
    cmap = np.ones(seq_len)
    v.add_colormap(values=[3], style='energy')
    #v.add_aux_BP(1, 10, color='red')
    v.savefig(out_fig)
    #v.show()


    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #df = pd.read_csv("/Users/timshel/NanoporeAnalysis/DashML/ShapeMap/varna_RNAse_P.csv")
    #df['Position'] = df['Position'] + 1
    #varnaapi.load_config('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    # sequence = 'GUUAAUCAUGCUCGGGUAAUCGCUGCGGCCGGUUUCGGCCGUAGAGGAAAGUCCAUGCUCGCACGGUGCUGAGAUGCCCGUAGUGUUCGUGCCUAGCGAAUCCAUAAGCUAGGGCAGCCUGGCUUCGGCUGGGCUGACGGCGGGGAAAGAACCUACGUCCGGCUGGGAUAUGGUUCGAUUACCCUGAAAGUGCCACAGUGACGGAGCUCUAAGGGAAACCUUAGAGGUGGAACGCGGUAAACCCCACGAGCGAGAAACCCAAAUGAUGGUAGGGGCACCUUCCCGAAGGAAAUGAACGGAGGGAAGGACAGGCGGCGCAUGCAGCCUGUAGAUAGAUGAUUACCGCCGGAGUACGAGGCGCAAAGCCGCUUGCAGUACGAAGGUACAGAACAUGGCUUAUAGAGCAUGAUUAACGUC'
    # ss = "(((((((((((((((((((((.(((((((((....)))))))))...[.[[.[[[[[(((((((((.(((...).)))))).....((((((((((((........)))))))((((((((((....))))))))).)((.(((((((((((((..((((.....))))).))))))).))))).....(((((............((((((((....)))))))).........)))..))))))))))))))...((((.........))))...((((((((((((.(...)))))))))))))(((((((........)))))))........)))))))((((((((((..(((.....))).......))))))))))......]]]]]]]].).)))))))))))))..."
    # v = varnaapi.Structure(structure=ss, sequence=sequence)
    #v = varnaapi.FileDraw('/Users/timshel/NanoporeAnalysis/DashML/VARNA/RNAseP_as-in-original-paper.bpseq')
    #v.savefig("RNAse_Pp.png")


##### Auto Generate Figure of Native Secondary Structure or Dominant Cluster
#### highlight conserved regions on figure
### manual generation for main article images using bpseq
def get_varna_putative_plots(df, cons_ss, cons_bp, metric, cluster=0):
    # use control structure if exists, else dominant cluster
    ss = df['control_secondary'].unique()[0]
    if (ss is None) or (len(ss) <= 0):
        df = df.loc[df['cluster'] == cluster]
        if len(df) <= 0:
            return
        ss = df['secondary'].unique()[0]
    lid = df['LID'].unique()[0]
    seq = df['contig'].unique()[0]
    seq_len = len(df['sequence'].unique()[0])
    sequence = df['sequence'].unique()[0]
    get_vplot(lid, seq, seq_len, sequence, ss, cons_ss, cons_bp, metric)
    return



####
### cons_threshold based on dominant cluster percentage
### region size default 3
### km_threshold based on SHAPE thresholds .4
def get_conserved_regions(lids, region_size=3, gap=1, centroid_threshold=None):
    ### get centroids and secondary structures
    df = dbsel.select_centroidz(lids)
    #TODO save regions for analysis, label with region #s
    # region size can be option
    df_clust = dbsel.select_max_clusters(lids)

    #### hamming
    df_hamming = df.loc[df['method'] == 'hamming']
    #dominant cluster
    dt = df_clust.loc[df_clust['method'] == 'hamming']
    clust_num = dt['cluster_size'].sum()
    # average of top 5 clusters as threshold
    std = dt.nlargest(5, 'cluster_size')['cluster_size'].std()/clust_num
    max_clust_percent = dt.nlargest(5, 'cluster_size')['cluster_size'].mean()/clust_num - np.abs(std/2)
    dt = dt.max()
    max_clust = dt['cluster']
    cons_ss, cons_bp = get_conserved_hamming(df_hamming,cons_thresh=(max_clust_percent if centroid_threshold==None
                                    else centroid_threshold), region_size=region_size, gap=gap)
    # # plot native or dominant secondary structure with conserved regions highlighted in varna
    get_varna_putative_plots(df_hamming, cons_ss, cons_bp, metric='hamming', cluster=max_clust)


    #### kmeans
    df_km = df.loc[df['method']=='kmeans']
    # dominant cluster
    dt = df_clust.loc[df_clust['method'] == 'kmeans']
    clust_num = dt['cluster_size'].sum()
    std = dt.nlargest(5, 'cluster_size')['cluster_size'].std() / clust_num
    max_clust_percent = (dt.nlargest(5, 'cluster_size')['cluster_size'].mean() / clust_num) - np.abs(std*1.5)
    dt = dt.max()
    max_clust = dt['cluster']
    cons_ss, cons_bp = get_conserved_kmeans(df_km,cons_thresh=(max_clust_percent if centroid_threshold==None
                                else centroid_threshold),region_size=region_size, gap=gap, km_threshold=.4)
    # # plot native or dominant secondary structure with conserved regions highlighted in varna
    get_varna_putative_plots(df_km, cons_ss, cons_bp, metric='kmeans', cluster=max_clust)
    return


#lids = '37'
#get_conserved_regions(lids)
