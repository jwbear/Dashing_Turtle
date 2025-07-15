import sys
import os.path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerPathCollection
from sklearn.neighbors import LocalOutlierFactor
import DashML.Database_fx.Insert_DB as dbins


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#### source ####
dmso_sequences = None
acim_sequences = None
seq = None


def get_metric(dfr):
    #### aggregate counts of Predict, read_depth,
    #### another decider based on the percentage of modified reads ###
    dfr = dfr[['position', 'contig', 'Predict']]
    dfr['Predict'] = dfr['Predict'].astype('category')
    dfr = dfr.groupby(['position', 'contig', 'Predict'], observed=False).size().unstack(fill_value=0)
    dfr.reset_index(inplace=True)
    dfr['read_depth'] = dfr[-1] + dfr[1]
    dfr['percent_modified'] = dfr[-1] / dfr['read_depth']
    #### predict modification based on percent modified read depth ####
    mean = dfr['percent_modified'].mean()
    dfr['Predict'] = np.where(dfr['percent_modified'] > mean, -1, 1)
    return dfr
def update_legend_marker_size(handle, orig):
    "Customize size of the legend marker"
    handle.update_from(orig)
    handle.set_sizes([20])


def get_novelty():
    dmso = dmso_sequences
    acim = acim_sequences

    def novelty_signal():
        ### train on all dmso ####
        dftrain = dmso[['event_level_mean']]
        X_train = dftrain.to_numpy()
        dfx = acim[['event_level_mean']]
        X = dfx.to_numpy()

        # fit the model for novelty detection (novelty=True)
        clf = LocalOutlierFactor(n_neighbors=5, novelty=True, contamination=.4)
        clf.fit(X_train)
        # DO NOT use predict, decision_function and score_samples on X_train as this
        # would give wrong results but only on new unseen data (not used in X_train),
        # e.g. X_test, X_outliers or the meshgrid
        y_pred_test = clf.predict(X)
        # n_errors = (y_pred_test != ground_truth).sum()
        X_scores = clf.negative_outlier_factor_

        #df = pd.DataFrame({"Shape_Map": ground_truth, "Lof_Novelty": y_pred_test, "BaseType":np.array(dft["BaseType"])})
        acim['predict_signal'] = y_pred_test
        acim['varna_signal'] = np.where(y_pred_test == -1, 1, 0)
        #dft.to_csv(save_path + seq + "_lof_signal.csv")
        #local metric
        # dft = get_metric(dft)
        # sp = save_path + seq + "_lof_signal_metric.csv"
        # mx.get_Metric(dft, seq, sp)

    def novelty_dwell():
        ### train on all dmso ####
        dftrain = dmso[['event_length']]
        X_train = dftrain.to_numpy()
        dfx = acim[['event_length']]
        X = dfx.to_numpy()

        # fit the model for novelty detection (novelty=True)
        clf = LocalOutlierFactor(n_neighbors=5, novelty=True, contamination=.4)
        clf.fit(X_train)
        # DO NOT use predict, decision_function and score_samples on X_train as this
        # would give wrong results but only on new unseen data (not used in X_train),
        # e.g. X_test, X_outliers or the meshgrid
        y_pred_test = clf.predict(X)

        # df = pd.DataFrame({"Shape_Map": ground_truth, "Lof_Novelty": y_pred_test, "BaseType":np.array(dft["BaseType"])})
        #dft = acim[['position', 'contig', 'read_index']]
        acim['predict_dwell'] = y_pred_test
        acim['varna_dwell'] = np.where(y_pred_test == -1, 1, 0)


    novelty_signal()
    novelty_dwell()
    acim = acim[['LID', 'contig', 'read_index', 'position','predict_signal', 'predict_dwell', 'varna_dwell', 'varna_signal']]
    dbins.insert_lof(acim)
