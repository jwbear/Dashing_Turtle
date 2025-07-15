import sys
import platform
import re
import numpy as np
import pandas as pd
from DashML.Predict import Predict_BPP as bpp
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel


### Note: Predict > 3 is noisier with more coverage, better for average predictions
### Predict > 4 stabilizes in the deconvolution offers improved predictions
#TODO: get all non reactive signle molecule base pair probabilities
# and get intermolecular non reactive bp probabilities for caomparison

#### TODO ## could use dmso reactivities ???, need to calculate dmso reactivities and integrate

#### Paper GMM positional induced clusters reflect predictions, nice
# TODO: optimize average and indiv reads, both should use maximal clusters for calculation
### TODO: predict based on cluster size

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

MAX_THREADS = None
bpp.MAX_THREADS = MAX_THREADS

def scale_reactivities(reactivities):
    min = reactivities.min()
    max = reactivities.max()
    smin = 0
    smax = 2

    reactivities = ((reactivities - min) / (max - min)) * (smax - smin) + smin
    return reactivities

def get_mods(lids, continue_reads=False, vienna=False):
    # get sequence data from db
    df = dbsel.select_predict(lids)
    df.fillna(value=0, inplace=True)
    df = df.groupby(['LID', 'contig', 'read_index', 'position', 'Sequence'], as_index=False)[[ 'Predict_BC',
       'Predict_Signal', 'Predict_Dwell', 'Predict_Lofd', 'Predict_Lofs',
       'Predict_Gmm']].agg(lambda x: x.mode().iloc[0]).reset_index(drop=True)

    if len(df) <= 0:
        raise Exception("No predictions found. Please run predictions.")
        return


    # -1 is outlier, 1 is inlier
    df['Predict_BC'] = np.where(df['Predict_BC'] == -1, 1, 0)
    df['Predict_Signal'] = np.where(df['Predict_Signal'] == -1, 1, 0)
    df['Predict_Dwell'] = np.where(df['Predict_Dwell'] == -1, 1, 0)
    df.loc[(df['Predict_Dwell'] == 1) & (df['Sequence'] != 'G'), "Predict_Dwell"] = 3
    df.loc[(df['Predict_Dwell'] == 1) & (df['Sequence'] == 'G'), "Predict_Dwell"] = 2
    df.loc[df['Predict_Dwell'] == 0, "Predict_Dwell"] = -2
    df['Predict_Lofd'] = np.where(df['Predict_Lofd'] == -1, 1, 0)
    df['Predict_Lofs'] = np.where(df['Predict_Lofs'] == -1, 1, 0)
    df['Predict_Gmm'] = np.where(df['Predict_Gmm'] == -1, 1, 0) #1,1,0?
    df.fillna(value=0, inplace=True)
    #print(df.head())
    print(len(df))

    # calculates average reactivity over all reads
    # predict is single score, above threshold
    # no probabilities
    # should use maximal clusters for this scoring metric
    # Least optimal method
    def mean_read(df):
        df['Reactivity'] = df['Predict_BC'] + df['Predict_Signal'] + df['Predict_Dwell'] + \
                            df['Predict_Lofs'] + df['Predict_Lofd'] + df['Predict_Gmm']
        #dfr = df[['position', 'contig', 'read_index', 'Sequence', 'Reactivity']]
        df['VARNA'] = np.where(df['Reactivity'] >=6, 1, 0)
        df['Predict'] = np.where(df['VARNA'] == 1, -1, 1)
        #mx.get_Metric(df, seq + "_mean_")
        # TODO add probabilities, opt
        dbins.insert_rx_full(df)

    #mean_read(df)

    def read_depth(df):
        print("Read_Depth")
        #### aggregate counts of Predict, read_depth,
        #### another decider based on the percentage of modified reads ###

        ##### Prediction per read, all reads
        df['Reactivity'] = df['Predict_BC'] + df['Predict_Signal'] + df['Predict_Dwell'] + \
                           df['Predict_Lofs'] + df['Predict_Lofd'] + df['Predict_Gmm']
        df['VARNA'] = np.where(df['Reactivity'] >= 6, 1, 0)
        df['Predict'] = np.where(df['VARNA'] == 1, -1, 1)
        # percent reactivities for shape in full set
        df['Reactivity_score'] = df['Reactivity']/8
        #adjust read predict based on probabilities
        if vienna==True:
            if continue_reads==False:
                dbins.insert_read_depth_full_clear(df)
                dfx = bpp.get_predict_probabilities(df, [lids], .95, reactivity=True,
                                                    continue_reads=False)
            else:
                dfx = bpp.get_predict_probabilities(df, [lids], .95, reactivity=True,
                                              continue_reads=True)
        else:
            if continue_reads==False:
                #delete
                dbins.insert_read_depth_full_clear(df)
                #insert
                dbins.insert_read_depth_full(df)
            else:
                dbins.insert_read_depth_full(df)


        ##### Prediction percent modified and base pairing probability, averaged
        dfr = df[['LID','position', 'contig', 'read_index', 'Reactivity']]
        df_rd = dbsel.select_readdepth(lids)
        dfr['Predict'] = df['Predict'].astype('category')
        dfr = (dfr.groupby(['LID', 'position', 'contig', 'Predict'], observed=False).size().
               unstack(fill_value=0).reset_index())
        dfr = dfr.merge(df_rd, on=['LID', 'position'], how='left')
        if -1 not in dfr.columns: #no modifications
            dfr[-1]=0
        dfr['percent_modified'] = dfr[-1]/dfr['read_depth']
        # input for RNAfold shape reactivity
        dfr['RNAFold_Shape_Reactivity'] = scale_reactivities(dfr['percent_modified'])
        mean = dfr['percent_modified'].mean()
        std = dfr['percent_modified'].std()
        #print("mean ", mean)


        #### Predict modification based on percent modified read depth ####
        #### Averaged over all transcripts based on read depth
        ### Most accurate correlates with cluster findings!!!!
        #TODO: certainty = read depth/ tot reads, variance = unmod prediction/tot umod reads
        #TODO gmm by read and check values
        dfr['Predict'] = np.where(dfr['percent_modified'] > mean + std/3, -1, 1)
        ### adjust prediction with base pairing probabilities
        dfr.fillna(value=0, inplace=True)
        dfr = bpp.get_predict_probabilities(dfr, [lids], .95, reactivity=True,
                                            continue_reads=continue_reads)
        dfr.rename(columns={-1: 'Out_num', 1: 'In_num'}, inplace=True)
        dfr['VARNA'] = np.where(dfr['Predict'] == -1, 1, 0)
        #print(dfr.head())
        #insert to db
        dbins.insert_read_depth(dfr)
        return

    read_depth(df)

# putative sequences designated in db
# sequences = ['RNAse_P', "cen_3'utr", "cen_3'utr_complex", 'cen_FL', 'cen_FL_complex',
#              "ik2_3'utr_complex", 'ik2_FL_complex', 'T_thermophila', 'ik2_FL', 'HCV',
#              "ik2_3'utr", 'HSP70_HSPA1A']

##### test
#get_mods(('43'),continue_reads=False)
#sys.exit(0)
