import sys, re
import platform
import os.path
import traceback
import numpy as np
import pandas as pd
import math
import scipy.signal
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn import metrics
from sklearn.metrics import RocCurveDisplay
from sklearn.linear_model import LinearRegression
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
from sklearn.cluster import DBSCAN
from sklearn.neighbors import KernelDensity
import DashML.Database_fx.Insert_DB as dbins


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


#### source ####
dmso_sequences = None
acim_sequences = None
contig = None


def unit_vector_norm(x):
    # unit vector normalization
    x = np.where(x != 0, ((x- np.min(x))/ (np.max(x)- np.min(x))), 0)
    #x = stats.zscore(x)
    return x

def ksm(df):
    bc_acim = df
    bc_dmso = dmso_sequences

    c = ['basecall_reactivity']
    xr = bc_dmso[c].to_numpy().flatten()  # unit_vector_norm(bc_dmso[c].to_numpy()).flatten()
    yr = bc_acim[c].to_numpy().flatten()  # unit_vector_norm(bc_acim[c].to_numpy()).flatten()
    end = len(xr)
    ##### manual cdf to measure ks #####
    # CDF
    # CDF(x) = "number of samples <= x"/"number of samples"
    x1 = np.sort(xr)
    y1 = np.sort(yr)

    def ecdf(x, v):
        res = np.searchsorted(x, v, side='right') / x.size
        return res

    delta_ecdf = np.subtract(yr, xr)


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

    p = find_peaks(delta_ecdf, plateau_size=[0, 10])

    #real values of delta_ecdf normalized
    peaks_value = unit_vector_norm(delta_ecdf)

    #set peaks in vector of all bases
    s1 = np.ones(len(delta_ecdf))
    s1[p[0]] = -1
    # s2 = np.zeros(len(delta_ecdf))
    # for i,n in enumerate(p[1].get('peak_heights')):
    #     if i in p[0]:
    #         s2[i] = n

    #print(sequence)
    dfc = pd.DataFrame()
    #print(p[1].get('peak_heights'))
    #dfc['peak_weight'] = s2
    dfc['LID'] = df['LID']
    dfc['position'] = df['position']
    dfc['contig'] = contig
    dfc['is_peak'] = s1
    dfc['peak_height'] = delta_ecdf
    dfc['insertion'] = bc_dmso['insertion'].astype('float32').to_numpy() - bc_acim['insertion'].astype('float32').to_numpy()
    dfc['mismatch'] = bc_dmso['mismatch'].astype('float32').to_numpy() - bc_acim['mismatch'].astype('float32').to_numpy()
    dfc['deletion'] = bc_dmso['deletion'].astype('float32').to_numpy() - bc_acim['deletion'].astype('float32').to_numpy()
    dfc['quality'] = bc_dmso['quality'].astype('float32').to_numpy() - bc_acim['quality'].astype('float32').to_numpy()
    dfc['basecall_reactivity'] = bc_dmso['basecall_reactivity'].astype('float32').to_numpy() - bc_acim['basecall_reactivity'].astype('float32').to_numpy()
    dfc['aligned_reads'] = bc_dmso['aligned_reads'].astype('float32').to_numpy() - bc_acim['aligned_reads'].astype('float32').to_numpy()
    #dfc['is_peak2'] = np.where(dfc['peak_height'] > .008, 1, 0)
    #dfc.to_csv(sequence + "_weightcompare.csv")
    #print(dfc.head(100))
    return dfc

def get_bc_reactivity_peaks():
    # basecall error file
    acim = acim_sequences
    #print(acim.columns)

    dft = pd.DataFrame()
    for lid in acim['LID'].unique():
        df = acim[acim['LID']==lid]
        peak_weight = ksm(df=df)
        if peak_weight is not None:
            dft = pd.concat([dft, peak_weight], ignore_index=True)


    dft['Predict'] = np.where(dft['is_peak'] == -1, -1, 1)
    dft['VARNA'] = np.where(dft['is_peak'] == -1, 1, 0)
    dbins.insert_basecall(dft)
    return
