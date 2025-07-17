import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc

# true_positive_rate = [0.4, 0.8,  0.8, 1.0, 1]
# false_positive_rate = [0, 0,  0, 0.2, 0.8]
#
#
#### SHAPE MAP
# def aroc(file):
#     df = pd.read_csv(file)
#     df.drop(df[df.Predict == 0].index, inplace=True)
#     #print(df.columns)
#     y_true = np.where(df["base_type"]=='B', 1, -1).flatten()
#     y_score = df["Predict"].to_numpy().flatten()
#     df = pd.DataFrame(columns=['x', 'y'])
#     df['x'] = y_score
#     df['y'] = y_true
#     tp = len(df[(df['x']==-1) & (df['y']==-1)])
#     fp = len(df[(df['x'] == -1) & (df['y'] == 1)])
#     tn = len(df[(df['x'] == 1) & (df['y'] == 1)])
#     fn = len(df[(df['x'] == 1) & (df['y'] == -1)])
#     print(tp, fp, tn, fn)
#     ppv = tp / (tp + fp)
#     acc = (tp + tn) / (tp + tn + fp + fn)
#     print('ppv ', ppv)
#     print('acc ', acc)
#     tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
#     roc_auc = auc(fpr, tpr)
#     #print(roc_auc)
#     return tpr, fpr, thresholds, roc_auc
#
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/HCV_map.csv')
# plt.plot(x, y, label=f'HCV (AUC = {ra:.2f})', color='orange')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/T_thermophila_map.csv')
# plt.plot(x, y, label=f'T_thermophila (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/RNase_P_map.csv')
# plt.plot(x, y, label=f'RNase_P (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/E_coli_map.csv')
# plt.plot(x, y, label=f'E_coli (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/hc16_map.csv')
# plt.plot(x, y, label=f'hc16_ligase (AUC = {ra:.2f})')
# #plt.plot(fpr, tpr, 'o-', label='DT', color='orange')
# plt.plot([0, 1], [0, 1], '--', color='grey')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('AUC-ROC Curve SHAPE-MaP Data')
# plt.legend()
# plt.savefig(fname="shapemap.png", dpi=600)
# #plt.show)
#
#
# plt.clf()
#
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/HCV_ce.csv')
# plt.plot(x, y, label=f'HCV (AUC = {ra:.2f})', color='orange')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/T_thermophila_ce.csv')
# plt.plot(x, y, label=f'T_thermophila (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/RNAseP_ce.csv')
# plt.plot(x, y, label=f'RNase_P (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/E_coli_ce.csv')
# plt.plot(x, y, label=f'E_coli (AUC = {ra:.2f})')
# x,y, z, ra = aroc('/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/hc16_ce.csv')
# plt.plot(x, y, label=f'hc16_ligase (AUC = {ra:.2f})')
# #plt.plot(fpr, tpr, 'o-', label='DT', color='orange')
# plt.plot([0, 1], [0, 1], '--', color='grey')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('AUC-ROC Curve SHAPE-CE Data')
# plt.legend()
# plt.savefig(fname="shapece.png", dpi=600)
# #plt.show)
#
# plt.clf()
def stats(x, y):
    df = pd.DataFrame(columns=['x', 'y'])
    df['x'] = x
    df['y'] = y
    tp = len(df[(df['x']==-1) & (df['y']==-1)])
    fp = len(df[(df['x'] == -1) & (df['y'] == 1)])
    tn = len(df[(df['x'] == 1) & (df['y'] == 1)])
    fn = len(df[(df['x'] == 1) & (df['y'] == -1)])
    print(tp, fp, tn, fn)
    ppv = tp / (tp + fp)
    acc = (tp + tn) / (tp + tn + fp + fn)
    print('ppv ', ppv)
    print('acc ', acc)

def hcvroc(file, file1):
    df = pd.read_csv(file)
    df.drop(df[df.position > 376].index, inplace=True)
    #df.drop(df[df.position < 12].index, inplace=True)
    y_true = np.where(df["metric"] == 'B', 1, -1).flatten()
    df = pd.read_csv(file1)
    df['position'] = df['position'] - 1
    #df.drop(df[df.position < 12].index, inplace=True)
    y_score = df["Predict"].to_numpy().flatten()
    stats(y_score, y_true)
    tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return tpr, fpr, thresholds, roc_auc

def rnaseproc(file, file1):
    df = pd.read_csv(file)
    df.drop(df[df.position > 412].index, inplace=True)
    #df.drop(df[df.position < 12].index, inplace=True)
    y_true = np.where(df["metric"] == 'B', 1, -1).flatten()
    df = pd.read_csv(file1)
    df['position'] = df['position'] - 1
    #df.drop(df[df.position < 12].index, inplace=True)
    y_score = df["Predict"].to_numpy().flatten()
    stats(y_score, y_true)
    tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return tpr, fpr, thresholds, roc_auc

def tthermophilaroc(file, file1):
    df = pd.read_csv(file)
    df.drop(df[df.position > 381].index, inplace=True)
    #df.drop(df[df.position < 12].index, inplace=True)
    y_true = np.where(df["metric"] == 'B', 1, -1).flatten()
    df = pd.read_csv(file1)
    df['position'] = df['position'] - 1
    #df.drop(df[df.position > 387].index, inplace=True)
    y_score = df["Predict"].to_numpy().flatten()
    stats(y_score, y_true)
    tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return tpr, fpr, thresholds, roc_auc

def ecoliroc(file, file1):
    df = pd.read_csv(file)
    #df.drop(df[df.position > 362].index, inplace=True)
    #df.drop(df[df.position < 12].index, inplace=True)
    y_true = np.where(df["metric"] == 'B', 1, -1).flatten()
    df = pd.read_csv(file1)
    df['position'] = df['position'] - 1
    df.drop(df[df.position > 361].index, inplace=True)
    y_score = df["Predict"].to_numpy().flatten()
    stats(y_score, y_true)
    tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return tpr, fpr, thresholds, roc_auc

def hc16roc(file, file1):
    df = pd.read_csv(file)
    #df.drop(df[df.position > 381].index, inplace=True)
    #df.drop(df[df.position < 12].index, inplace=True)
    y_true = np.where(df["metric"] == 'B', 1, -1).flatten()
    df = pd.read_csv(file1)
    df['position'] = df['position'] - 1
    df.drop(df[df.position > 336].index, inplace=True)
    y_score = df["Predict"].to_numpy().flatten()
    stats(y_score, y_true)
    tpr, fpr, thresholds = roc_curve(y_true, y_score, pos_label=-1)
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return tpr, fpr, thresholds, roc_auc

x,y, z, ra = hcvroc('/DashML/Landscape/Article_Images/HCV_metric.csv',
                   '/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/HCV_read_depth.csv')
plt.plot(x, y, label=f'HCV (AUC = {ra:.2f})', color='orange')
x,y, z, ra = tthermophilaroc('/DashML/Landscape/Article_Images/T_thermophila_metric.csv',
                             '/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/T_thermophila_read_depth.csv')
plt.plot(x, y, label=f'T_thermophila (AUC = {ra:.2f})')
x,y, z, ra = rnaseproc('/DashML/Landscape/Article_Images/Rnasep_metric.csv',
                  '/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/RNAse_P_read_depth.csv')
plt.plot(x, y, label=f'RNase_P (AUC = {ra:.2f})')
x,y, z, ra = ecoliroc('/DashML/Landscape/Article_Images/Ecoli_metric.csv',
                  '/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/Ecoli_read_depth.csv')
plt.plot(x, y, label=f'E_coli (AUC = {ra:.2f})')
x,y, z, ra = hc16roc('/DashML/Landscape/Article_Images/hc16_metric.csv',
                  '/Users/timshel/structure_landscapes/DashML/Deconvolution/Landscape/Article_Images/hc16_read_depth.csv')
plt.plot(x, y, label=f'hc16_ligase (AUC = {ra:.2f})')
plt.plot([0, 1], [0, 1], '--', color='grey')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('AUC-ROC Curve DT Data')
plt.legend()
plt.savefig(fname="dt.png", dpi=600)
#plt.show)


### ppv
