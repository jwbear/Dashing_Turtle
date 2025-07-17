import sys, re
import os.path
import traceback
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
import DashML.Database_fx.Insert_DB as dbins


#### source ####
acim_sequences = None
contig = None


def get_metric(dfr):
    #### aggregate counts of Predict, read_depth,
    #### another decider based on the percentage of modified reads ###
    dfr = dfr[['LID', 'position', 'contig', 'Predict']]
    dfr['Predict'] = dfr['Predict'].astype('category')
    dfr = dfr.groupby(['LID', 'position', 'contig', 'Predict'], observed=False).size().unstack(fill_value=0)
    dfr.reset_index(inplace=True)
    dfr['read_depth'] = dfr[-1] + dfr[1]
    dfr['percent_modified'] = dfr[-1] / dfr['read_depth']
    #### predict modification based on percent modified read depth ####
    mean = dfr['percent_modified'].mean()
    dfr['Predict'] = np.where(dfr['percent_modified'] > mean, -1, 1)
    print(dfr.columns)
    return dfr

def ksm(delta):
    delta_ecdf = delta.to_numpy().flatten()
    end = len(delta_ecdf)

    log_dens = []
    n = 0
    width = 2
    for i in range(0,len(delta_ecdf)):
        m = n
        n = n + width
        if n < len(delta_ecdf) and m < len(delta_ecdf):
            kde = KernelDensity(kernel="gaussian", bandwidth="silverman").fit(delta_ecdf.reshape(-1, 1)[m:n])
            l= kde.score_samples(delta_ecdf.reshape(-1,1)[m:n])
            log_dens.append(l)
        else:
            kde = KernelDensity(kernel="gaussian", bandwidth="silverman").fit(delta_ecdf.reshape(-1, 1)[m-width:])
            l = kde.score_samples(delta_ecdf.reshape(-1, 1)[m:])
            log_dens.append(l)
            break;

    tmp = []
    for l in log_dens:
        for v in l:
            tmp.append(v)
    #print(tmp)
    #kde peaks
    #delta_ecdf = tmp

    #sklearn signal processing peaks
    #delta_ecdf = dfx.unit_vector_norm(delta_ecdf)
    #print(delta_ecdf)
    #sys.exit(0)

    p = find_peaks(delta_ecdf, plateau_size=[0, 2])


    #set peaks in vector of all bases
    peaks = np.ones(len(delta_ecdf))
    peaks[p[0]] = -1
    #peaks[p[1]['left_edges']] = -1
    #peaks[p[1]['right_edges']] = -1
    # print(p[1]['left_edges'])
    # print(f"{np.abs(tmp[60:80])} {shape_outliers[60:80]} \n{s1[60:80]}")
    #sys.exit()

    pred_outliers = np.where((peaks == -1), True, False).sum()
    #print("Total Predicted Outlier Sites: ", pred_outliers)

    delta_ecdf = delta_ecdf
    return peaks, delta_ecdf


def get_reactivity_peaks():
    df = acim_sequences
    #print(df.columns)

    peaks, height = ksm(df['delta_dwell'])
    # print(np.mean(height))
    if peaks is not None:
        df['predict_dwell'] = peaks
        df['varna_dwell'] = np.where(peaks == -1, 1, 0)

    peaks, height = ksm(df['delta_signal'])
    # print(np.mean(height))
    if peaks is not None:
        df['predict_signal'] = peaks
        df['varna_signal'] = np.where(peaks == -1, 1, 0)

    print("Insert peaks..........")
    dbins.insert_peaks(df)

    print("Dwell Metric........")
    dt = df
    dt.rename(columns={'delta_dwell': 'delta', 'predict_dwell': 'Predict', 'varna_dwell': 'VARNA'}, inplace=True)
    dt = dt[['LID', 'read_index', 'contig', 'position', 'delta', 'Predict', 'VARNA']]


    print("Signal Metric........")
    ds = df
    ds.rename(columns={'delta_signal': 'delta', 'predict_signal': 'Predict', 'varna_signal': 'VARNA'}, inplace=True)
    ds = ds[['LID', 'read_index', 'contig', 'position', 'delta', 'Predict', 'VARNA']]
    return
