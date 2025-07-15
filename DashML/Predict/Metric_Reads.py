import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import warnings
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_Metrics_DB as dbsel
from concurrent.futures import ThreadPoolExecutor

MAX_THREADS = 10000
def get_statistics(x, y):
    tn = np.sum(np.where((x==1) & (y==1), 1, 0))
    tp = np.sum(np.where((x==-1) & (y==-1), 1, 0))
    fp = np.sum(np.where((x==-1) & (y==1), 1, 0))
    fn = np.sum(np.where((x==1) & (y==-1), 1, 0))
    accuracy = (tp + tn ) / (tp + tn + fp + fn)
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    ppv = tp / (tp + fp)
    p1, pearsons= pearsonr(x, y)
    m1, mannwhit= wilcoxon(x, y)
    aggreement = np.sum(np.where(x==y, 1, 0))/len(x)
    # print(aggreement)
    # print(tp, tn, fp, fn)
    # print(accuracy, sensitivity, specificity, ppv)
    # print(pearsons)
    # print(mannwhit)
    return tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement


### loops/ bulges predicted
def get_secondary_metric(secondary_labels, algorithm_labels):
    # count loops/bulges
    loop = 0
    loop_detect = 0
    bulge = 0
    bulge_detect = 0
    loop_coverage = 0
    bulge_coverage = 0
    loop_tot = 0
    bulge_tot = 0

    b = secondary_labels
    alg = algorithm_labels
    i = 0
    lenb = len(b)
    while i < lenb:
        is_detect = False
        if b[i] == 'Loop':
            loop = loop + 1
            while b[i] == 'Loop':
                if alg[i] == -1:
                    is_detect = True
                    loop_coverage = loop_coverage + 1
                loop_tot = loop_tot + 1
                i = i + 1
                if i >= lenb:
                    break
            if is_detect:
                loop_detect = loop_detect + 1

        elif b[i] == 'Bulge':
            bulge = bulge + 1
            while b[i] == 'Bulge':
                if alg[i] == -1:
                    is_detect = True
                    bulge_coverage = bulge_coverage + 1
                bulge_tot = bulge_tot + 1
                i = i + 1
                if i >= lenb:
                    break
            if is_detect:
                bulge_detect = bulge_detect + 1
        else:
            if i < len(b):
                i = i + 1
            else:
                b = None

    #print("Loops", loop, loop_detect, loop_coverage, loop_tot)
    #print("Bulges", bulge, bulge_detect, bulge_coverage, bulge_tot)
    return loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot

def control_metric(df, algorithm, threshold, metric, vienna=False):
    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    dfx = dft.groupby(['LID', 'read_index']).size().reset_index()
    dfx = dfx[['LID', 'read_index']]
    ids = dfx.loc[:, 'LID':'read_index'].values.tolist()
    #print('len ' + str(len(ids)))

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn', 'pearsons', 'mannwhit',
                               'aggreement'])

    def calc_metric(lid, read_index):
        dfx = dft.loc[(dft['LID']==lid) & (dft['read_index']==read_index)]
        if vienna:
            x = np.where(dfx['predict'] >= threshold, 1, -1)
        else:
            x = dfx['predict'].to_numpy()
        y = dfx['control'].to_numpy()
        tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
            get_statistics(x=x, y=y))
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]
        return

    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        for lid, read in ids:
            executor.submit(calc_metric, lid, read)
    return dt

def structure_metric(df, ycol, algorithm, threshold, metric, vienna=False):
    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn', 'pearsons', 'mannwhit',
                               'aggreement'])
    dfx = df.groupby(['LID', 'read_index']).size().reset_index()
    dfx = dfx[['LID', 'read_index']]
    ids = dfx.loc[:, 'LID':'read_index'].values.tolist()

    def calc_metric(lid, read_index):
        dfx = df.loc[(df['LID'] == lid) & (df['read_index'] == read_index)]
        if vienna:
            x = np.where(dfx['predict'] >= threshold, 1, -1)
        else:
            x = dfx['predict'].to_numpy()
        if ycol == 'base_type':
            y = np.where(dfx['base_type'] == 'S', -1, 1)
        elif ycol == 'metric':
            y = np.where(dfx['metric'] == 'B', 1, -1)
        tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
            get_statistics(x=x, y=y))
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        for lid, read in ids:
            executor.submit(calc_metric, lid, read)
    return dt

def secondary_metric(df, algorithm, threshold, metric, vienna=False):
        dt = pd.DataFrame(
            columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                     'bases_covered', 'bases_total'])
        dfx = df.groupby(['LID', 'read_index']).size().reset_index()
        dfx = dfx[['LID', 'read_index']]
        ids = dfx.loc[:, 'LID':'read_index'].values.tolist()

        def calc_metric(lid, read_index):
            dfx = df.loc[(df['LID'] == lid) & (df['read_index'] == read_index)]
            x = dfx['metric'].to_numpy()
            if vienna:
                y = np.where(dfx['predict'] >= threshold, 1, -1)
            else:
                y = dfx['predict'].to_numpy()
            ## compare to control structure secondary features loops/bulges
            loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
                get_secondary_metric(secondary_labels = x, algorithm_labels=y))
            metric = 'loops'
            dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, loop, loop_detect, loop_coverage,
                               loop_tot]
            metric = 'bulges'
            dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, bulge, bulge_detect, bulge_coverage,
                               bulge_tot]

        with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
            for lid, read in ids:
                executor.submit(calc_metric, lid, read)

        return dt

#### Algorithm Results ######

#### GMM ####
# TODO OPT Extend GMM Analysis to a per read version for comparison
def gmm(lids):
    print("GMM Averaged, Read Index -1")

#### Basecall Peaks ####
# TODO OPT Extend BCto a per read version for comparison
def basecall(lids):
    print("BC Averaged, Read Index -1")

#### Signal Peaks ####

def signal_peaks(lids):
    algorithm = "signal peaks"
    threshold = 2 ##peak width

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_signal_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Signal Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df,algorithm, threshold, metric)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_signal_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df,algorithm, threshold, metric)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_signal_control(lids, ['map','ce'])
    dt = pd.concat([dt, control_metric(df,algorithm, threshold, metric)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_signal_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt,structure_metric(df,'metric', algorithm, threshold, metric)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))

#### dwell Peaks ####
def dwell_peaks(lids):
    algorithm = "dwell peaks"
    threshold = 2 ##peak width

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_dwell_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Dwell Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_dwell_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_dwell_control(lids, ['map', 'ce'])
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_dwell_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt, structure_metric(df, 'metric', algorithm, threshold, metric)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))


#### Lof Signal ####
def lof_signal(lids):
    algorithm = "lof signal"
    threshold = 0 ##parameters only no threshold

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_lofs_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: LOFS Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_lofs_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_lofs_control(lids, ['map', 'ce'])
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_lofs_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt, structure_metric(df, 'metric', algorithm, threshold, metric)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))


#### Lof Dwell ####
def lof_dwell(lids):
    algorithm = "lof dwell"
    threshold = 0 ##parameters only no threshold

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_lofd_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: LOFD Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_lofd_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_lofd_control(lids, ['map', 'ce'])
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_lofd_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt, structure_metric(df, 'metric', algorithm, threshold, metric)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))


#### Vienna Comparison #####
def vienna(lids):
    algorithm = "vienna"
    threshold = .95 ##default parameters

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn', 'pearsons', 'mannwhit',
                               'aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_vienna_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Vienna Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric, vienna=True)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_vienna_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric, vienna=True)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_vienna_control(lids, ['map', 'ce'])
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric, vienna=True)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_vienna_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric, vienna=True)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt, structure_metric(df, 'metric', algorithm, threshold, metric, vienna=True)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric, vienna=True)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))


#### Reactivity Results ######

#### Reactivity ####
def read_depth_full(lids):
    algorithm = "read depth full"
    threshold = 6 ## of 8 total point metrics

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_rdf_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: RDF Predictions not calculated for sequence.')
        return
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    metric = 'ce'
    df = dbsel.select_rdf_control(lids, 'ce')
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    # xor
    metric = 'shape_xor'
    df = dbsel.select_rdf_control(lids, ['map', 'ce'])
    dt = pd.concat([dt, control_metric(df, algorithm, threshold, metric)], ignore_index=True)

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_rdf_structure(lids)
    dt = pd.concat([dt, structure_metric(df, 'base_type', algorithm, threshold, metric)], ignore_index=True)

    metric = 'structure acim'
    dt = pd.concat([dt, structure_metric(df, 'metric', algorithm, threshold, metric)], ignore_index=True)

    # insert dt
    dbins.insert_metric(dt)

    # count loops/bulges
    dt = secondary_metric(df, algorithm, threshold, metric)

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt.head(10))

def read_depth(lids):
    print("Read Depth Averaged, Read Index -1")


#lids = 36

def get_metrics(lids):
    ## averaged not app
    #gmm([lids])
    #basecall([lids])
    #read_depth([lids])
    ## by read index
    signal_peaks([lids])
    dwell_peaks([lids])
    lof_signal([lids])
    lof_dwell([lids])
    read_depth_full([lids])
    vienna([lids])
    return
