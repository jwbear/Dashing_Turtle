import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc


df = pd.DataFrame({'Sequence':['HCV IRES', 'T. thermophila', 'RNASE P', 'E. coli tmRNA', 'hc16 ligase'],
                   'Shape-Map':[78, 47, 68,74,81],'DT':[85,62,72,70,78], 'Shape-CE' :[75,40,76,55,70],  })
df.index = df['Sequence']
ax = df.plot.bar(rot=0, color=['blue', 'purple', 'grey'])
ax.set_ylabel('Pos. Pred. Value %')
plt.title('Positive Predictive Value')
plt.savefig(fname="ppv.png", dpi=600)
#plt.show)
