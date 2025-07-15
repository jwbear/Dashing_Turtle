import sys, os
import platform
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import DashML.Database_fx.Select_DB as dbsel
import DashML.Database_fx.Insert_DB as dbins

# lids = '50'
# df = dbsel.select_read_depth_ave(lids)
# df.rename(columns={'Rnafold_shape_reactivity':'Reactivity'}, inplace=True)
# seq_name = 'AverageReactivity_' + df['contig'].unique()[0] + '_' + str.replace(lids, ',', '-')
# ax = df.plot.bar(x='position', y='Reactivity',xlabel='Position', ylabel='Reactivity', rot=90, figsize=(12, 6))
# # Set the y-axis minimum
# ax.set_ylim(0,2)
# ax.set_xlim(5,110)
# ax.set_xticks(np.arange(5, 110, 5))
# ax.set_title('Average Predicted reactivity')
# figname = os.path.join('Figures/' + seq_name + '.png')
# plt.savefig(figname)
# #plt.show)
#
#
# df.rename(columns={'Base_pair_prob':'Base Pairing Prob.'}, inplace=True)
# seq_name = 'AverageBasePairingProbability_' + df['contig'].unique()[0] + '_' + str.replace(lids, ',', '-')
# ax = df.plot.bar(x='position', y='Base Pairing Prob.',xlabel='Position', ylabel='Probability', rot=90, figsize=(12, 6))
# # Set the y-axis minimum
# ax.set_ylim(0,1)
# ax.set_xlim(5,110)
# ax.set_xticks(np.arange(5, 110, 5))
# ax.set_title('Maximum Predicted Base Pairing Probabilities (ViennaRNA')
# figname = os.path.join('Figures/' + seq_name + '.png')
# plt.savefig(figname)
# #plt.show)

dtr = pd.DataFrame.from_dict({
                "contig": ['tester'],
                "sequence": ['acgtu'],
                "secondary": ['(...)'],
                "experiment": ['nmr'],
                "sequence_name": ['polly'],
                "sequence_len": [5],
                "temp": [37],
                "is_modified": [0],
                "type1": ['dmso'],
                "type2": ['acim'],
                "complex": [0],  # TODO: Add complex input if needed
                "run": [1]
            }, orient='columns')
print(dtr.head())
lid = dbins.insert_library(dtr)
print(lid)
