import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import platform
import pandas as pd
from platformdirs import user_documents_path
from pathlib import Path
from DashML.Predict import Peak_Analysis_BC as bc
from DashML.Predict import Peak_Analysis as peak
from DashML.Predict import Gmm_Analysis2 as gmx
from DashML.Predict import Predicts as predict_mods
from DashML.Predict import Predict_BPP as predict_probs
from DashML.Predict import Lof
from DashML.Database_fx import Select_DB as dbsel


print(user_documents_path())
output = user_documents_path() / "DTLandscape_Output"
output.mkdir(parents=True, exist_ok=True)

#### Paper GMM positional induced clusters reflect predictions, nice
# TODO: add to Predict, + for indicate
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.options.mode.chained_assignment = None

def get_lids(seq, temp, type1, type2, complex):
    ### Source Data ####
    # TODO: Add modified modified case, generalizing when type1==type2
    # generate new lid, duplicated Mod table for that contig with new id in stored proc
    # take into account multiple runs require diff ids
    # then run sub predictions normally
    if type1==type2:
        ### predicts on unmodified data, same data sets
        acim_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type2, complex=complex)
        dmso_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type1, complex=complex)
    else:
        acim_sequences = dbsel.select_mod(seq, temp, type1=type1, type2=type2, complex=complex)
        dmso_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type1)
    unmod_seq_ids = dmso_sequences['LID'].unique()
    unmod_seq_ids = str(unmod_seq_ids).replace("[", "").replace("]", "")
    mod_seq_ids = acim_sequences['LID'].unique()
    mod_seq_ids = str(mod_seq_ids).replace("[", "").replace("]", "")

    return mod_seq_ids
####Predict by Desc #####
def predict(seq, temp, type1, type2, complex):

    ### Source Data ####
    # TODO: Add modified modified case, generalizing when type1==type2
    # generate new lid, duplicated Mod table for that contig with new id in stored proc
    # take into account multiple runs require diff ids
    # then run sub predictions normally

    if type1==type2:
        ### predicts on unmodified data, same data sets
        acim_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type2, complex=complex)
        dmso_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type1, complex=complex)
    else:
        acim_sequences = dbsel.select_mod(seq, temp, type1=type1, type2=type2, complex=complex)
        dmso_sequences = dbsel.select_unmod(seq, temp, type1=type1, type2=type1)
    unmod_seq_ids = dmso_sequences['LID'].unique()
    unmod_seq_ids = str(unmod_seq_ids).replace("[", "").replace("]", "")
    mod_seq_ids = acim_sequences['LID'].unique()
    mod_seq_ids = str(mod_seq_ids).replace("[", "").replace("]", "")

    if (len(acim_sequences) <= 0) or (len(dmso_sequences) <= 0):
        raise Exception("No sequences found. Please load sequence data.")
        return

    ### TODO: remove extra parameters from subtasks, and update subtask inserts

    #### Basecall Analysis
    def bc_analysis():
        print("Running BC Anaylysis for " + seq)
        ## bc pre_run diff source data
        bc.acim_sequences = dbsel.select_bc_mod(mod_seq_ids)
        bc.dmso_sequences = dbsel.select_bc_unmod(unmod_seq_ids)
        if (len(bc.acim_sequences) <= 0) or (len(bc.dmso_sequences) <= 0):
            raise Exception("Basecall Data Insufficient. Please load basecall data.")
            return
        bc.contig = seq
        bc.get_bc_reactivity_peaks()
    bc_analysis()


    #### peak analysis
    def peak_analysis():
        print("Running Peak Anaylysis for " + seq)
        if type1==type2:
            # TODO fix for mod mod
            peak.acim_sequences = dbsel.select_peaks_unmod(mod_seq_ids, unmod_seq_ids)
        else:
            peak.acim_sequences = dbsel.select_peaks(mod_seq_ids, unmod_seq_ids)
        peak.contig = seq
        peak.get_reactivity_peaks()
    peak_analysis()


    ##### lof
    def lof_analysis():
        print("Running Lof Anaylysis for " + seq)
        Lof.acim_sequences = acim_sequences
        Lof.dmso_sequences = dmso_sequences
        Lof.seq = seq
        Lof.get_novelty()
    lof_analysis()


    ##### gmm
    def gmm_analysis():
        print("Running GMM Anaylysis for " + seq)
        gmx.acim_sequences = acim_sequences
        gmx.dmso_sequences = dmso_sequences
        gmx.seq = seq
        gmx.seq_ids = mod_seq_ids
        gmx.positional_gmm()
    gmm_analysis()

    return mod_seq_ids


####Predict by LID #####
def predict_lids(unmod_lids, mod_lids):

    ### Source Data ####

    ### predicts on unmodified data, same data sets
    acim_sequences = dbsel.select_mod_lid(mod_lids)
    dmso_sequences = dbsel.select_unmod_lid(unmod_lids)


    if (len(acim_sequences) <= 0) or (len(dmso_sequences) <= 0):
        seq_err = ("No sequences found. Please load sequence data. modified {} unmodified {}".
                   format(len(acim_sequences), len(dmso_sequences)))
        raise Exception (seq_err)
        print("Error:", seq_err, file=sys.stderr)
        return None

    unmod_seq_ids = dmso_sequences['LID'].unique()
    unmod_seq_ids = str(unmod_seq_ids).replace("[", "").replace("]", "")
    mod_seq_ids = acim_sequences['LID'].unique()
    mod_seq_ids = str(mod_seq_ids).replace("[", "").replace("]", "")
    seq = str(acim_sequences['contig'].unique()[0]).replace("[", "").replace("]", "")

    ### TODO: remove extra parameters from subtasks, and update subtask inserts

    #### Basecall Analysis
    def bc_analysis():
        print("Running BC Anaylysis for " + seq)
        ## bc pre_run diff source data
        bc.acim_sequences = dbsel.select_bc_mod(mod_seq_ids)
        bc.dmso_sequences = dbsel.select_bc_unmod(unmod_seq_ids)
        if (len(bc.acim_sequences) <= 0) or (len(bc.dmso_sequences) <= 0):
            sys.error("Basecall Data Insufficient. Please load basecall data.")
            sys.exit(0)
        bc.contig = seq
        bc.get_bc_reactivity_peaks()
    bc_analysis()


    #### peak analysis
    def peak_analysis():
        print("Running Peak Anaylysis for " + seq)
        peak.acim_sequences = dbsel.select_peaks(mod_seq_ids, unmod_seq_ids)
        peak.contig = seq
        peak.get_reactivity_peaks()
    peak_analysis()


    ##### lof
    def lof_analysis():
        print("Running Lof Anaylysis for " + seq)
        Lof.acim_sequences = acim_sequences
        Lof.dmso_sequences = dmso_sequences
        Lof.seq = seq
        Lof.get_novelty()
    lof_analysis()


    ##### gmm
    def gmm_analysis():
        print("Running GMM Anaylysis for " + seq)
        gmx.acim_sequences = acim_sequences
        gmx.dmso_sequences = dmso_sequences
        gmx.seq = seq
        gmx.seq_ids = mod_seq_ids
        gmx.positional_gmm()
    gmm_analysis()

    return mod_seq_ids

#### Predict Reactivity for All Reads and Averaged #####
def get_mods(lids, continue_reads, vienna=False):
    predict_mods.MAX_THREADS = 100  # vary by platform speed up computation
    predict_mods.get_mods(lids, continue_reads=continue_reads, vienna=vienna)


#### Get Vienna Probabilities for Exisitng Predicitons #####
def base_prob(lids, threshold=.95, continue_reads=False, gui=None):
    predict_probs.MAX_THREADS = 100  # vary by platform speed up computation
    predict_probs.get_predict_probability(lids, threshold=threshold, continue_reads=continue_reads)


#### Get Graphs of Predictions for GUI #####
def get_graph_ave(lids):
    df = dbsel.select_read_depth_ave(lids)
    df.rename(columns={'Rnafold_shape_reactivity': 'Reactivity'}, inplace=True)
    seq_name = 'AverageReactivity_' + df['contig'].unique()[0] + '_' + str.replace(str(lids), ',', '-')
    filtered_df = df[(df['position'] >= 0) & (df['position'] <= 100)]
    ax = filtered_df.plot.bar(x='position', y='Reactivity', legend=False, figsize=(8, 2), width=0.8)
    # Set the y-axis minimum
    ax.set_ylim(0, 2)
    # Set ticks only every 5th position
    positions = df['position'].values
    tick_indices = np.arange(0, len(positions), 5)  # every 5th bar
    ax.set_xticks(tick_indices)
    ax.set_xticklabels([positions[i] for i in tick_indices], rotation=45, fontsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.set_xlabel('Position', fontsize=8)
    ax.set_ylabel('Reactivity', fontsize=8)
    ax.set_title('Average Predicted Reactivity', fontsize=9)
    subdir = output / "Figures"
    subdir.mkdir(parents=True, exist_ok=True)
    figname1 = subdir / f"{seq_name}.png"
    figname1 = str(figname1.resolve())
    fig = ax.get_figure()
    fig.tight_layout()
    plt.savefig(figname1)
    plt.close(fig)
    ##plt.show)

    # print current data
    subdir = output / "Data"
    subdir.mkdir(parents=True, exist_ok=True)
    f_name = ('AvePredictedReactivity_' +
              df['contig'].unique()[0] + '_' + str.replace(str(lids), ',', '-')) + '.csv'
    file_path = subdir / f_name
    file_path = str(file_path.resolve())
    df[['LID', 'contig', 'position', 'Reactivity']].to_csv(file_path, index=False)

    # graph of base pairing probabilities
    df = df.rename(columns={'Base_pair_prob': 'Base Pairing Prob.'})
    filtered_df = filtered_df.rename(columns={'Base_pair_prob': 'Base Pairing Prob.'})
    seq_name = 'AveMaxBasePairingProbability_' + df['contig'].unique()[0] + '_' + str.replace(str(lids), ',', '-')
    ax = filtered_df.plot.bar(
        x='position', y='Base Pairing Prob.', legend=False,
        figsize=(8, 2), width=0.8
    )
    # Set the y-axis minimum
    ax.set_ylim(0, 1)
    # Set ticks only every 5th position
    positions = df['position'].values
    tick_indices = np.arange(0, len(positions), 5)  # every 5th bar
    ax.set_xticks(tick_indices)
    ax.set_xticklabels([positions[i] for i in tick_indices], rotation=45, fontsize=8)
    # Set font size on y-axis ticks
    ax.tick_params(axis='y', labelsize=8)
    ax.set_xlabel('Position', fontsize=8)
    ax.set_ylabel('Probability', fontsize=8)
    ax.set_title('Ave. Max Base Pairing Prob. (ViennaRNA)', fontsize=9)
    subdir = output / "Figures"
    subdir.mkdir(parents=True, exist_ok=True)
    figname2 = subdir / f"{seq_name}.png"
    figname2 = str(figname2.resolve())
    fig = ax.get_figure()
    fig.tight_layout()
    plt.savefig(figname2)
    plt.close(fig)
    # #plt.show)

    # print current data
    subdir = output / "Data"
    subdir.mkdir(parents=True, exist_ok=True)
    f_name = ('AveMaxBasePairingProbability_' +
              df['contig'].unique()[0] + '_' + str.replace(str(lids), ',', '-')) + '.csv'
    file_path = subdir / f_name
    file_path = str(file_path.resolve())
    #print(df.columns)
    df[['LID', 'contig', 'position', 'Base Pairing Prob.']].to_csv(file_path, index=False)

    return figname1, figname2




#### Run Unmodified & Modified Data Separately
#### Different Runs of same contig are grouped

### USER DEFINED: LID SELECTED BY DROP DOWN OF POPOPULATED BY DB ####
#seq = 'HCV' #contig in FAFSA/nanopolish
#temp = 37 #37, 42 other
#type1 = 'dmso' #unmodified condition, drop down
#type2 = 'acim' #modified condition, drop down
#complex = 0 #is contig in complex
threshold = .95
#unmod_lids = 52
#lids = 37


### Calculate modifications
# run if done separetly (default)
def run_predict(unmod_lids, lids, continue_reads=False, vienna=False):
    try:
        ### Get LIDs for modification parameters
        #mod_lids = get_lids(seq=seq, temp=temp, type1=type1, type2=type2, complex=complex)

        #### Run Sub Predictions
        #predict(seq=seq, temp=temp, type1=type1, type2=type1, complex=complex)
        predict_lids(unmod_lids, lids)

        ### Run Combined Predictions by Read and Averaged
        # checkboxes:
        #   vienna determines if bpp should be calculated at same time, long running, default is false
        #   continue_reads allows process to be restarted, default is True
        #   False deletes any previous data and restarts bpp
        get_mods(lids, continue_reads=continue_reads, vienna=vienna)

        ### Graphs and Prints Averaged Results ####
        image_path1, image_path2 = get_graph_ave(lids)

        ### METRICS
        ### get metrics on averaged predictions ###
        #metric_avg.get_metrics(mod_lids)
        ### get metrics per read ###
        # TODO add import

        return image_path1, image_path2
    except Exception as err:
        raise Exception("Prediction error: " + str(err))
        return None, None
    finally:
        print("Predictions Complete")


#run_predict(26,43, continue_reads=False, vienna=True)

### Calculate modification probabilities over all predictions Vienna (optional)
# run if done separetly (default)
#base_prob(lids, threshold, continue_reads=False)
