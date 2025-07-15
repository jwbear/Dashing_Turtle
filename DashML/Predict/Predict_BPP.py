import sys, os
import traceback
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from DashML.Predict import Predict_Fold as pbpp
from concurrent.futures import ThreadPoolExecutor
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

### get probabilities for a single sequence for dmso, eg no cluster reactivities
### could use dmso reactivities ???, need to calculate dmso reactivities

### Combines probabilities extracted from rnacofold and get max bpp for each base
### Used in Predict
# get max probabilities for single molecule folding of A and AA interactions

MAX_THREADS = 100

def scale_reactivities(reactivities):
    min = reactivities.min()
    max = reactivities.max()
    smin = 0
    smax = 2

    reactivities = ((reactivities - min) / (max - min)) * (smax - smin) + smin
    return reactivities

### TODO incomplete
def get_prob_pool_sep(df, lids, threshold, continue_reads=False):
    # get probabilities for each lid-lid interactions
    ids = []
    ssids = []

    if continue_reads== True:
        dfr = dbsel.select_continued_reads(lids)
        dfx = df.groupby(['LID', 'read_index']).size().reset_index()
        dfx = dfx.merge(dfr, on=['LID', 'read_index'], how='left')
        dfx['completed'].fillna(0, inplace=True)
        dfx = dfx.loc[dfx['completed'] == 0]
        dfx = dfx[['LID', 'read_index']]
        ids = dfx.loc[:, 'LID':'read_index'].values.tolist()
        del dfr
        del dfx
    else:
        dfx = df.groupby(['LID', 'read_index']).size().reset_index()
        dfx = dfx[['LID', 'read_index']]
        ids = dfx.loc[:, 'LID':'read_index'].values.tolist()
        del dfx


    def get_rx(lid, read):
        try:
            dt = df.loc[(df['LID'] == int(lid)) & (df['read_index'] == int(read))]
            reactivity = dt['Reactivity_score'].to_numpy().flatten()
            reactivity = scale_reactivities(reactivity)
            ssid = pbpp.get_probabilities(lid, lid, reactivity=reactivity, read=read)
        except Exception as e:
            raise Exception(str(e))
            print("rx exception",e)
            return



    def get_probpool():
        with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
            for lid, read in ids:
                executor.submit(get_rx, lid, read)

    get_probpool()

    # add bpp to rdf
    dbins.insert_read_depth_full_update(lids, threshold, rx_threshold=.7)
    return df

def get_prob_pool(df, lids, threshold, continue_reads=False):
    # get probabilities for each lid-lid interactions
    ids = []
    ssids = []
    df_ps = pd.DataFrame()


    if continue_reads== True:
        dfr = dbsel.select_continued_reads(lids)
        dfx = df.groupby(['LID', 'read_index']).size().reset_index()
        dfx = dfx.merge(dfr, on=['LID', 'read_index'], how='left')
        dfx['completed'].fillna(0, inplace=True)
        dfx = dfx.loc[dfx['completed'] == 0]
        dfx = dfx[['LID', 'read_index']]
        ids = dfx.loc[:, 'LID':'read_index'].values.tolist()
        del dfr
        del dfx
    else:
        dfx = df.groupby(['LID', 'read_index']).size().reset_index()
        dfx = dfx[['LID', 'read_index']]
        ids = dfx.loc[:, 'LID':'read_index'].values.tolist()
        del dfx


    def get_rx(lid, read):
        dt = df.loc[(df['LID'] == int(lid)) & (df['read_index'] == int(read))]
        reactivity = dt['Reactivity_score'].to_numpy().flatten()
        reactivity = scale_reactivities(reactivity)
        ssid = pbpp.get_probabilities(lid, lid, reactivity=reactivity, read=read)
        df_ps = dbsel.select_max_structure_probabilities(ssid)
        if len(df_ps) > 0:
            df_ps['position'] = df_ps['position'] - 1
            if 'base_pair_prob' in dt.columns:
                dt = dt.drop(columns='base_pair_prob')
            dt = dt.merge(df_ps, on=['LID', 'position', 'read_index'], how='left')
            dt['completed'] = 1
            ### Adjust Predict Based on Threshold ####
            dt['Predict'] = np.where(
                ((dt['Predict'] == -1) & (dt['Reactivity_score'] < .7) & (dt['base_pair_prob'] > threshold)), 1,
                dt['Predict'])
        if 'base_pair_prob' in dt.columns:
            dt['base_pair_prob'].fillna(0, inplace=True)
        try:
            dbins.insert_read_depth_full(dt)
        except Exception as e:
            raise Exception({e})
            return
        return

    def get_probpool():
        with ThreadPoolExecutor(max_workers=100) as executor:
            for lid, read in ids:
                executor.submit(get_rx, lid, read)

    get_probpool()

    return df

def get_predict_probabilities(df, lids, threshold, reactivity=False, continue_reads=False):
    df.sort_values(by=['LID', 'position'], inplace=True)
    # get probabilities for each lid-lid interactions
    df_ps = pd.DataFrame()

    if reactivity == True:
        #read_depth_full
        if 'Reactivity_score' in df.columns: #reactivity for each read read_depth_full
           df = get_prob_pool(df, lids, threshold, continue_reads=continue_reads)

        #read_depth
        elif 'RNAFold_Shape_Reactivity' in df.columns: # averaged reactivity read_depth, single value
            for lid in lids:
                reactivity = df.loc[df['LID']==int(lid), 'RNAFold_Shape_Reactivity'].to_numpy().flatten()
                ssid = pbpp.get_probabilities(lid, lid,  reactivity=reactivity, read=-1)
                df_ps = pd.concat([df_ps, dbsel.select_max_structure_probabilities(ssid)])
            df_ps['position'] = df_ps['position'] - 1
            df['read_index'] = -1
            df = df.merge(df_ps, on=['LID', 'position', 'read_index'], how='left')

            ### Adjust Predict Based on Threshold ####
            df['Predict'] = np.where(((df['Predict'] == -1) & (df['RNAFold_Shape_Reactivity']<.8)
                                      & (df['base_pair_prob'] > threshold)), 1, df['Predict'])
            df['base_pair_prob'].fillna(0, inplace=True)

    else: # no reactivity for unmodified and averaged if unmodified reactivity not calculated
        for lid in lids:
            ssid = pbpp.get_probabilities(lid, lid, reactivity=None, read=-2)
            df_ps = pd.concat([df_ps, dbsel.select_max_structure_probabilities(ssid)])
        df_ps['position'] = df_ps['position'] -1
        df['read_index'] = -2
        df = df.merge(df_ps, on=['LID', 'position', 'read_index'], how='left')
        print(df)

        ### Adjust Predict Based on Threshold ####
        df['Predict'] = np.where(
            ((df['Predict'] == -1) & (df['Reactivity_score'] < .7) & (df['base_pair_prob'] > threshold)), 1,
            df['Predict'])
        df['base_pair_prob'].fillna(0, inplace=True)

    return df

#separate run for deriving vienna probabilities for predicted lids
def get_predict_probability(lids, threshold=.95, continue_reads=False):
    try:
        if isinstance(lids, int):
            lids = str(lids)
            lids = [lids]

        df = dbsel.select_read_depth_full_all(lids)
        df.sort_values(by=['LID', 'position'], inplace=True)
        # get probabilities for each lid-lid interactions
        df_ps = pd.DataFrame()

        # check if reactivities exist in read depth full tables, if not return
        if 'Reactivity_score' in df.columns:  # reactivity for each read read_depth_full
            df = get_prob_pool_sep(df, lids, threshold, continue_reads=continue_reads)
        else: # no reactivity for unmodified and averaged if unmodified reactivity not calculated
            for lid in lids:
                ssid = pbpp.get_probabilities(lid, lid, reactivity=None, read=-2)
                if ssid is not None:
                    df_ps = pd.concat([df_ps, dbsel.select_max_structure_probabilities(ssid)])
                    df_ps['position'] = df_ps['position'] -1
            df['read_index'] = -2
            if len(df_ps) > 0:
                df = df.merge(df_ps, on=['LID', 'position', 'read_index'], how='left')
            print(df)

            ### Adjust Predict Based on Threshold ####
            df['Predict'] = np.where(
                ((df['Predict'] == -1) & (df['Reactivity_score'] < .7) & (df['base_pair_prob'] > threshold)), 1,
                df['Predict'])
            df['base_pair_prob'].fillna(0, inplace=True)
    except Exception as e:
        raise Exception(str(e))
        print("Probability Error ", e)
        return
    finally:
        return
