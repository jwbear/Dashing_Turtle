import sys
import numpy as np
import pandas as pd
import warnings
import math
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_Metrics_DB as dbsel

def get_statistics(x, y):
    try:
        (tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons,
         mannwhit, aggreement) = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
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
        print(aggreement)
        print(tp, tn, fp, fn)
        print(accuracy, sensitivity, specificity, ppv)
        print(pearsons)
        print(mannwhit)
        return tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement
    except Exception as e:
        sys.stderr.write(str(e))
    finally:
        vz = ['tn', 'tp', 'fp', 'fn', 'accuracy', 'sensitivity', 'specificity', 'ppv', 'pearsons', 'mannwhit', 'aggreement']
        for v in vz:
            print(vars()[v])
            if math.isnan(vars()[v]):
                vars()[v] = 0

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

    print(len(secondary_labels), len(algorithm_labels))
    if (len(secondary_labels)<=1) or (len(algorithm_labels)<=1):
        return loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot

    b = secondary_labels
    alg = algorithm_labels
    i = 0
    while i < len(b)-1:
        is_detect = False
        if b[i] == 'Loop':
            loop = loop + 1
            while b[i] == 'Loop':
                if alg[i] == -1:
                    is_detect = True
                    loop_coverage = loop_coverage + 1
                loop_tot = loop_tot + 1
                i = i + 1
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
            if is_detect:
                bulge_detect = bulge_detect + 1
        else:
            if i < len(b)-1:
                i = i + 1
            else:
                b = None

    print("Loops", loop, loop_detect, loop_coverage, loop_tot)
    print("Bulges", bulge, bulge_detect, bulge_coverage, bulge_tot)
    return loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot

#### BASElINE CONTROL RESULTS ######
#structure compare shape,ce, xor to structure
#### get metrics for control experiments
def shape_map(lids):
    algorithm = "shape_map"
    read_index = -1
    threshold = .4

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### get control data, change types/df if additional controls added
    df = dbsel.select_shape_control(lids, "map", 'ce')

    #### compare to ce
    metric = 'ce'
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    #### count loops/bulges
    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=x))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

def shape_ce(lids):
    algorithm = "shape_ce"
    read_index = -1
    threshold = .4

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### get control data, change types/df if additional controls added
    df = dbsel.select_shape_control(lids, "ce", 'map')

    #### compare to ce
    metric = 'map'
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    #### count loops/bulges
    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=x))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

def vienna(lids):
    algorithm = "vienna"
    read_index = -2 #indicated predicted structure without reactivity
    threshold = 0 #predicted mfe structure w/no reactivity

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn', 'pearsons', 'mannwhit',
                               'aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_gmm_control(lids, 'map')
    # remove 0 values, missing in control
    df = df.loc[df['gmm'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    metric = 'ce'
    df = dbsel.select_gmm_control(lids, 'ce')
    # remove 0 values, missing in control
    df = df.loc[df['control'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_gmm_control(lids, ['map', 'ce'])
    # remove 0 values, missing in control
    df = df.loc[df['control'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_gmm_structure(lids)
    x = df['gmm'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)
    print(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot = (
        get_secondary_metric(secondary_labels=df['metric'].to_numpy(), algorithm_labels=df['gmm'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, loop, loop_detect, loop_coverage, loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, bulge, bulge_detect, bulge_coverage,
                           bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)


#### Algorithm Results ######

#### GMM ####
def gmm(lids):
    algorithm = "gmm"
    read_index = -1
    threshold = .02 ###percent modified, relative modification rate of basecall reactivity

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_gmm_control(lids, 'map')
    # remove 0 values, missing in control
    df = df.loc[df['gmm'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_gmm_control(lids, 'ce')
    # remove 0 values, missing in control
    df = df.loc[df['control'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_gmm_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    df = df.loc[df['control'] != 0]
    x = df['gmm'].to_numpy()
    y = df['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_gmm_structure(lids)
    x = df['gmm'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)
    print(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['gmm'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### Basecall Peaks ####
def basecall(lids):
    algorithm = "basecall_peaks"
    read_index = -1
    threshold = 2 ##peak width

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_bc_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Basecall Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_bc_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_bc_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_bc_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### Signal Peaks ####
def signal_peaks(lids):
    algorithm = "signal peaks"
    read_index = -1
    threshold = 2 ##peak width

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_signal_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Signal Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_signal_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_signal_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_signal_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### dwell Peaks ####
def dwell_peaks(lids):
    algorithm = "dwell peaks"
    read_index = -1
    threshold = 2 ##peak width

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_dwell_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Dwell Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_dwell_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_dwell_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_dwell_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### Lof Signal ####
def lof_signal(lids):
    algorithm = "lof signal"
    read_index = -1
    threshold = 0 ##parameters only no threshold

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_lofs_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: LOFS Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_lofs_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_lofs_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_lofs_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### Lof Dwell ####
def lof_dwell(lids):
    algorithm = "lof dwell"
    read_index = -1
    threshold = 0 ##parameters only no threshold

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_lofd_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: LOFD Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_lofd_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_lofd_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_lofd_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

#### Vienna Comparison #####
def vienna(lids):
    algorithm = "vienna"
    read_index = -1
    threshold = .95 ##default parameters

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn', 'pearsons', 'mannwhit',
                               'aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_vienna_control(lids, 'map')
    if df['predict'].mean() <= 0:
        warnings.warn('Warning: ViennaRNA probabilities not calculated for sequence. '
                       'Run base pairing probabilities to continue...')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = np.where(dft['predict']>=threshold, 1, -1)
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    metric = 'ce'
    df = dbsel.select_vienna_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = np.where(dft['predict']>=threshold, 1, -1)
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_vienna_control(lids, ['map', 'ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = np.where(dft['predict']>=threshold, -1, 1)
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_vienna_structure(lids)
    x = np.where(df['predict']>=threshold, -1, 1)
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot = (
        get_secondary_metric(secondary_labels=df['metric'].to_numpy(),
                             algorithm_labels=np.where(df['predict']>=threshold, -1, 1)))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, loop, loop_detect, loop_coverage, loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, bulge, bulge_detect, bulge_coverage,
                           bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)


#### Reactivity Results ######

#### Reactivity ####
def read_depth_full(lids):
    algorithm = "read depth full"
    read_index = -1
    threshold = 6 ## of 8 total point metrics

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_rdf_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_rdf_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_rdf_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_rdf_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)

def read_depth(lids):
    algorithm = "read depth"
    read_index = -1
    threshold = .3 ## mean+std/3 of percent modified

    dt = pd.DataFrame(columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric','ppv', 'accuracy',
                               'sensitivity', 'specificity', 'tp', 'fp', 'tn', 'fn','pearsons', 'mannwhit','aggreement'])

    ### compare to experimental controls
    metric = 'shape_map'
    df = dbsel.select_rd_control(lids, 'map')
    if len(df) <= 0:
        warnings.warn('Warning: Predictions not calculated for sequence.')
        return

    # remove 0 values, missing in control
    dft = df.loc[df['predict'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index,threshold, metric, ppv, accuracy,sensitivity, specificity,
                  tp, fp,tn, fn, pearsons,mannwhit,aggreement]

    metric = 'ce'
    df = dbsel.select_rd_control(lids, 'ce')
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement =(
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # xor
    metric = 'shape_xor'
    df = dbsel.select_rd_control(lids, ['map','ce'])
    # remove 0 values, missing in control
    dft = df.loc[df['control'] != 0]
    x = dft['predict'].to_numpy()
    y = dft['control'].to_numpy()
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                tp, fp, tn, fn, pearsons, mannwhit, aggreement]


    ## compare to control structure
    metric = 'structure base type'
    df = dbsel.select_rd_structure(lids)
    x = df['predict'].to_numpy()
    y1 = np.where(df['base_type'] == 'S', -1, 1)
    y2 = np.where(df['metric'] == 'B', 1, -1)
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y1))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                  tp, fp, tn, fn, pearsons, mannwhit, aggreement]
    metric = 'structure acim'
    tn, tp, fp, fn, accuracy, sensitivity, specificity, ppv, pearsons, mannwhit, aggreement = (
        get_statistics(x=x, y=y2))
    for lid in lids:
        dt.loc[len(dt)] = [lid, algorithm, read_index, threshold, metric, ppv, accuracy, sensitivity, specificity,
                           tp, fp, tn, fn, pearsons, mannwhit, aggreement]

    # insert dt
    dbins.insert_metric(dt)

    dt = pd.DataFrame(
        columns=['LID', 'algorithm', 'read_index', 'threshold', 'metric', 'number_total', 'number_detected',
                 'bases_covered', 'bases_total'])

    ## compare to control structure secondary features loops/bulges
    loop, loop_detect, loop_coverage, loop_tot, bulge, bulge_detect, bulge_coverage, bulge_tot =  (
        get_secondary_metric(secondary_labels = df['metric'].to_numpy(), algorithm_labels=df['predict'].to_numpy()))

    metric = 'loops'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, loop, loop_detect, loop_coverage,loop_tot]

    metric = 'bulges'
    for lid in lids:
        dt.loc[len(dt)] = [lid,algorithm, read_index, threshold,metric, bulge, bulge_detect, bulge_coverage, bulge_tot]

    # insert dt
    dbins.insert_secondary_metric(dt)

    print(dt)
#mlid=36

def get_metrics(lids):
    gmm(lids)
    shape_map(lids)
    shape_ce(lids)
    basecall(lids)
    signal_peaks(lids)
    dwell_peaks(lids)
    lof_signal(lids)
    lof_dwell(lids)
    read_depth_full(lids)
    read_depth(lids)
    vienna(lids)
    return

#get_metrics([36])
