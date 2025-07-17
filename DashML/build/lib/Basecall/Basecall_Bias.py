import os.path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Modification_Bias:

    def __init__(self,reference_sequences, f_path, structure="", loop_structure="", dir_name="Default"):
        self.dir_name = dir_name
        self.reference_sequences = reference_sequences
        self.f_path = f_path
        #get base pairing data
        if structure:
            self.structure = structure
        #get loop structure data
        if 'single_loop_structure' in loop_structure:
            self.single_loop_structure = loop_structure.get('single_loop_structure')
        if 'paired_loop_structure' in loop_structure:
            self.paired_loop_structure = loop_structure.get('paired_loop_structure')
        ###### Create Output Folder ifne ##########
        self.save_path = self.f_path + 'Out_Structure/' + self.dir_name + '_Modification_Rates'
        print("Modification Output To: " + self.save_path)


    ########### Plot Modifications #################
    def parse_modifications(self, position, indel):
        df = pd.DataFrame()
        pymol = {}
        deletions = []
        insertions = []
        mismatch = []
        mismatches = []
        summary = []

        for key in position:
            sequence = [char for char in self.reference_sequences[key]]
            for i in range(len(position[key])):
                if type(indel[key][i]) is not int:
                    if indel[key][i][0] > 0:
                        deletions.append((position[key][i], indel[key][i][0]))
                    if indel[key][i][1] > 0:
                        insertions.append((position[key][i], indel[key][i][1]))
                    if indel[key][i][2] > 0:
                        mismatch.append((position[key][i], indel[key][i][2]))
                    if len(indel[key][i][3]) > 0:
                        #pos, number of mismatches,mismatch nucleotide,refseq nucleotide
                        mismatches.append((position[key][i], indel[key][i][3][0][1], indel[key][i][3][0][0],sequence[i]))
                    if (indel[key][i][0] > 0 or indel[key][i][1] > 0 or indel[key][i][2] > 0):
                        summary.append((position[key][i], indel[key][i][0], indel[key][i][1], indel[key][i][2]))
            ### if modifications then plot ########
            #if len(deletions) > 0:
                 #pltmod.plot_modification(deletions, seq_name=key, mod_type="Deletions", dir_name=dir_name)
            #if len(insertions) > 0:
                #pltmod.plot_modification(insertions, seq_name=key, mod_type="Insertions", dir_name=dir_name)
            #if len(mismatch) > 0:
                #pltmod.plot_modification(mismatch, seq_name=key, mod_type="Mismatch", dir_name=dir_name)
            if len(mismatches) > 0:
                 #### calculate mismatch bias ####
                 pos_data, df_tmp = self.mismatch_bias(mismatches, key)
                 df = df.append(df_tmp)
                 pymol[key] = self.pymol_parse(data=pos_data)
                 #pltmod.plot_mismatch(mismatches, seq_name=key, dir_name=dir_name)
            ##### Plot Summary of Modifications ############
            #if len(summary) > 0:
                #pltmod.plot_modification_summary(summary, dir_name=dir_name, seq_name=key)
            deletions = []
            insertions = []
            mismatch = []
            mismatches = []
            summary = []

        ######### Table of Modification Commands for Pymol ##########
        self.print_pymol_modifications(pymol)

        #### Table of Combined Modification Frequencies #####
        ####### Modifications by position ############
        df["Position"] = pd.to_numeric(df["Position"])
        df = df.reset_index(drop=True)
        #print(df.dtypes) #print(df.columns)
        df = df.loc[:, ["Sequence", "Position", "Mismatches","Number of Mismatches","Total Mismatches in Aligned Sequences", "Modified Bp Freq"]]
        df.sort_values(by=['Sequence', 'Position'], ascending=[True, True], inplace=True)
        self.print_df(df, f_suffix="_combined")
        ####### Modifications in sequence grouped by mismatch ##########
        df = df.groupby("Mismatches", as_index=False)["Number of Mismatches"].sum()
        df.insert(2, "Modified Bp Freq", (df.iloc[:,1]/df.iloc[:,1].sum()) * 100, allow_duplicates=True)
        #print(df)

    ######## parse pymol commmands for modified positions ##########
    def pymol_parse(self, data):
        # Prints list of modfications for pymol
        # TODO: auto generate pymol movie fx
        positions = ""
        for i in range(len(data)):
            positions += data[i] + "+"
        return positions

    ######## generate pymol commmands for modified positions ##########
    def print_pymol_modifications(self, pymol):
        # Prints list of modfications for pymol
        # TODO: auto generate pymol movie fx
        pymol_df = pd.DataFrame.from_dict(pymol, orient="index")
        pymol_df = pymol_df.reset_index()
        pymol_df.columns = ["Sequence", "Mismatch Positions"]
        self.print_df(pymol_df, f_suffix="_pymol")


    #### calculate mismatch bias per sequence ####
    def mismatch_bias(self, mismatches, curr_seq):
        base = re.compile('[A-Z]', re.IGNORECASE)
        bias = {}
        pymol = []

        data = np.array(mismatches)
        df = pd.DataFrame(mismatches)
        df.columns = ["Position", "Number of Mismatches", "Called Base", "Reference Base"]
        df["Position"] = pd.to_numeric(df["Position"])
        df['Called Base'] = df['Called Base'].str.replace('[^A-Z]', '', regex=True)
        df['Reference Base'] = df['Reference Base'].str.replace('[^A-Z]', '', regex=True)
        df['Mismatches'] = df['Reference Base'] + df['Called Base']
        df = df.drop(columns=['Called Base', 'Reference Base'])
        df.insert(3, "Unmodified Structure", df['Position'].apply(lambda x: 'S' if x in self.structure['single_strand_positions'] else 'BP'), allow_duplicates=True)

        #print(df)
        # if len(self.single_loop_structure) > 1:
        #     df_single_loop_structure = pd.DataFrame(self.single_loop_structure, columns=["Position", "Reference Base", "Loop Structure"])
        #     df_single_loop_structure = df_single_loop_structure.drop(["Reference Base"], axis =1)
        #     df = df.join(df_single_loop_structure.set_index('Position'), on='Position')
        #     df = df.fillna(0)
        #     #print(df_single_loop_structure)
        #     #print(df)
        df.sort_values(by=['Position', 'Mismatches', 'Unmodified Structure'], inplace=True)
        ######### Print Tables TODO: why??? ############
        self.print_df(df)
        ########## Plot Modification Data #############
        self.plot_df(df, curr_seq)

        self.print_df(df, curr_seq, "_Structural_Mismatches")

        #### total number of mismatches in aligned sequences * number of aligned reads
        # TODO: may want to include number of aligned reads
        tot = np.sum(data[:, 1].astype(np.int))
        #print("Tot: " + str(tot))
        #print(data[:,2])

        ######## compile string for pymol command #####
        self.pymol_parse(data[:, 0])

        #TODO: add position for strand to loop data[:, 0]
        for i in range(len(data[:, 3])):
            #remove artifacts from pysam indel notation
            ref_base = re.search(base, np.array2string(data[i, 3])).group()
            called_base = re.search(base, np.array2string(data[i, 2])).group()
            if ref_base != called_base:
                key = ref_base + called_base
                value = data[i,1].astype(np.int)
                if bias.get(key):
                    curr_val = bias.get(key)[0]
                else: curr_val = 0
                pos = data[i,0]
                if curr_val == None:
                    bias[key] = (value, pos)
                else:
                    bias[key] = (curr_val + value, pos)

        ######## Prepare Data Frame ##############
        df = pd.DataFrame.from_dict(bias, orient="index")
        df = df.reset_index()
        df.columns = ["Mismatches", "Number of Mismatches","Position"]
        df.insert(0, "Sequence", [curr_seq]*len(df), allow_duplicates=True)
        #print(df)
        df.insert(3, "Total Mismatches in Aligned Sequences", df[ "Number of Mismatches"].sum(), allow_duplicates=True)
        ### TODO: total mismatches at position percentage would be better
        df.insert(4, "Modified Bp Freq", df["Number of Mismatches"]/df["Total Mismatches in Aligned Sequences"] * 100, allow_duplicates=True)
        df = df.reset_index(drop=True)
        self.print_df(df, curr_seq, f_suffix="_modification_bias")
        return data[:, 0], df

        # TODO: add whether a modification occurs on single strand or on known base pairs, then later in/out helix
        #TODO: add function calculates bias vs guassian, vs unmodified tetrahymena
        #TODO: finally bias calculation for all modifications

    ######## Print Data Frame to CSV ##############
    def print_df(self,df, curr_seq="", f_suffix=""):

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        ###### output modification bias for file ###########
        if curr_seq:
            df.to_csv(self.save_path + "/" + curr_seq + "_" + self.dir_name + f_suffix + ".csv", index=False)
        else:
            df.to_csv(self.save_path + "/" + self.dir_name + f_suffix + ".csv", index=False)

    def plot_df(self, df, curr_seq):
        df['Number of Mismatches'] = (df['Number of Mismatches']/df['Number of Mismatches'].sum() ) * 100
        df = df.loc[df['Number of Mismatches'] > df['Number of Mismatches'].mean()]
        df['Mismatches'] = df['Mismatches'] + df['Unmodified Structure']
        #print(df)

        fig = df.plot(x='Position', y = 'Number of Mismatches', kind="bar", rot=45, title= curr_seq + " Mismatches by Structure", figsize=(12,8))
        # fig = df.plot(kind="bar", rot=90, title="Modification Rates by Position", figsize=(12,8), stacked=True)
        #### save plot #####
        i = 0
        score = df['Mismatches'].to_numpy()
        for p in fig.patches:
            fig.annotate(score[i], xy=(p.get_x(), p.get_height()))
            i += 1

        fig = plt.gcf()

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        figname = os.path.join(self.save_path + '/'+ curr_seq + '_' + 'Position_Structure_Modification_Rate' + '.png')
        fig.savefig(figname)
        #plt.show)
