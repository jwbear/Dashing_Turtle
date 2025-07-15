import sys, re
import os.path
import platform
import traceback
import numpy as np
import pandas as pd
import math
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import DashML.data_fx as dfx

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if platform.system() == 'Linux':
    ##### server #####
    data_path = "/home/jwbear/projects/def-jeromew/jwbear/Deconvolution/"
    save_path = "/home/jwbear/projects/def-jeromew/jwbear/Deconvolution/Out/"
else:
    data_path = sys.path[1] + "/DashML/Deconvolution/"
    save_path = sys.path[1] + "/DashML/Deconvolution/Latex/"


# get reactivities for native structures using hamming binary reactivities
# remove sections with no data for each cluster, only compare non-zero positions
def get_native_mfe(seq="HCV"):
    df_native = dfx.get_structure()
    df_native = df_native.loc[df_native['Sequence_Name']==seq]
    #print(df_native.head())
    df_native = df_native[['BaseType']]
    df_native['BaseType'] = np.where(df_native['BaseType']=='S', 1, 0)
    mm = df_native.to_numpy()
    #print(mm.shape)
    return mm.T

def latex_table(df, cols=None):
    cols = ['Sequence', 'Canonical Modification Rates', 'Non-Canonical Modification Rates']
    df = pd.DataFrame(columns=cols)
    latex = df.to_latex(index=False,float_format="{:.2f}".format)
    latex = str(latex).replace("_", "\\_")
    with open(save_path + "canonical.tex", 'w+') as f:
        f.write(latex)
