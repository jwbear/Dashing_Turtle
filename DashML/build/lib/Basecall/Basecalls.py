### Parse Alignnment File (BAM/SAM format) ######
#* Calculates dels, insertions, and mismatches
#* Calculates modification rates by base, and by sequence
#* Inputs: references sequences.fasta, sam file, bam file with baq calculated
#* Outputs: plots by sequences of dels, insertion, mismatches, incl. summary plot of all mods
#* normalized to the unit vector
#### NOTE pysam reads are referenced according to reference file usually in 5' - 3' direction

import sys
import os.path
import pysam
import statistics
import DashML.Basecall.Basecall_Plot as pltmod
from DashML.Basecall.Basecall_Bias import *
from collections import Counter



######### TODO: automate samtools commands ##############
######### TODO: automate alignment ##############
######### TODO: separate data output from algorithm ##########
# TODO: check for U/T conversions in mismatch T -> U should not be counted,
#  its not ref seq has T then mismatch wont count
# TODO: check all positions in signal analysis start with 0
# TODO: add minimap pairwise alignment with pdb structure sequence,
#  generate list of ins, del, start stop- maybe can do with PYSAM
#  pip install mappy- minimap2 library would like a list of contiguous matches
#   maybe add mismatch/ins/deletion notation in annotation file

class ParseAlignment:
    def __init__ (self, reference_sequences, f_path, f_reference, f_sam, f_bam, dir_name, save_path):
        ######### Set file paths and output names ##############
        self.reference_sequences = reference_sequences
        self.f_path = f_path
        self.f_reference = f_reference
        self.f_sam = f_sam
        self.f_bam = f_bam
        self.dir_name = dir_name
        self.save_path = save_path

        ### Create Pysam File Objects ###
        self.fafile = pysam.Fastafile(self.f_reference)
        self.samfile = pysam.AlignmentFile(self.f_sam, "r")
        self.bamfile = pysam.AlignmentFile(self.f_bam)

        ######### dictionary of summary modifications ##########
        self.sarr = {}  # total number of alignments by position
        self.baq = {}  # quality of base alignment
        self.consensus_pos = {}  # modification position in reference sequence, -1 if no reads
        self.indel = {}  # all modifications as tuple (ins, del, mismatch, mismatch base name)
        ######### Print Header ###############
        print("Processing following FASTA: " )
        print(self.bamfile.header)


    def get_base_modifications_by_read(self):
        bases = []
        try:
            for read_num, read in enumerate(self.bamfile.fetch()):
                #read id in sam files references to read_index in fastq
                sd = read.to_dict()
                #print("Read number: ", read_num)
                #print(sd.get('name'), sd.get('ref_name'))
                r = read.get_aligned_pairs(matches_only=False, with_seq=True)
                #print(r)
                start_read = read.query_alignment_start
                end_read = read.query_alignment_end
                #print("Read offset start-end: " + str(start_read) + '-' + str(end_read))
                #print("Read length:", end_read - start_read)
                # insertion would be represented by (read_pos, None)
                # deletion by (None, ref_pos)
                for index, i in enumerate(r):
                    # -2 no read
                    # 0 no change
                    # -1 deletion
                    # 1 insertion
                    # 2 mismatch
                    # read_num, position, error
                    if (i[0] == None) & (i[1] != None) & (i[2]!=None):
                        #print(str(i[1]) + ' Deletion')
                        bases.append([sd['ref_name'], read_num, i[1], -1])
                    elif (i[0] != None) & (i[1] == None) & (i[2]==None):
                        ref_offset = int(i[0])
                        if (ref_offset>=start_read) & (ref_offset<=end_read):
                            # check no position given, should be fine if consistent
                            # multiple insertions count as 1
                            #print(str(r[index-1][1]) + ' Insertion')
                            bases.append([sd['ref_name'], read_num, r[index-1][1], 1])
                    elif (i[0] != None) & (i[1] != None) & (i[2]!=None):
                        if i[2].islower():
                            #print(str(i[1]) + ' Mismatch')
                            bases.append([sd['ref_name'], read_num, i[1], 2])
                        else:
                            #print(str(i[1]) + ' ' + str(i[2]))
                            bases.append([sd['ref_name'], read_num, i[1], 0])
                #print(read.get_reference_sequence())
        except Exception as v:
            sys.stderr.write(str(v))
            sys.stderr.write("\nSequence failed\n")
            return

        # remove duplicates for insertion
        df = pd.DataFrame(bases, columns=['contig', 'read_index', 'position', 'basecall_error'])
        pad = []
        for r in df['read_index'].unique():
            dt = df.loc[df['read_index'] == r]
            min = dt['position'].min()
            max = dt['position'].max()
            seq = dt['contig'].unique()[0]
            seq_len = len(self.reference_sequences[seq])
            for p in np.arange(0, min):
                # pad front of sequence
                pad.append({'position': p, 'contig': seq, 'read_index': r, 'basecall_error':-2})
            for p in np.arange(max+1, seq_len):
                # pad end of sequence
                pad.append({'position': p, 'contig': seq, 'read_index': r, 'basecall_error': -2})
        df = pd.concat([df, pd.DataFrame(pad)], ignore_index=True)
        df = df.groupby(by=['contig', 'read_index', 'position'])['basecall_error'].max().reset_index()
        df = df.sort_values(by=['read_index', 'position'])
        ###### Create Output Folder ifne ##########
        save_path = self.save_path + self.dir_name + '_Modification_Reads/'
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        # save modifications by read
        df.to_csv( save_path + self.dir_name +'_modifications_by_read.csv', index=False)
        #check padding
        #print(df['position'].min(), df['position'].max())
        return

    def get_base_modifications(self):
        data = {} # dict of function values
        i = 0  #flag for loop
        ######## Create Pileups ##############
        iter = self.bamfile.pileup(fastafile=self.fafile)
        ########## Iterate Aligned Reads by Sequence (assumes MSA) #################
        while True:
            try:
                pileColumn = next(iter)
                curr_position = pileColumn.reference_pos
                #initialize first time only
                if i == 0:
                    curr_ref = pileColumn.reference_name
                    ########## temporary arrays store base information ##########
                    ref_len = self.bamfile.header.get_reference_length(curr_ref)
                    sna = [0] * ref_len  # number of alignments
                    bq = [0] * ref_len # quality of base alignment
                    pos = [-1] * ref_len # modification position in reference sequence
                    pos_indel = [-1] * ref_len # list of modifications position in aligned sequence
                    i+=1
                ##### iterate by reference name #############
                #### TODO: ref name not accesible in iterator??? ##############
                if curr_ref == pileColumn.reference_name:
                    ##### Number of Alignments #####
                    # num of alignments for position in qry, for all reads aligned to target
                    num_align = pileColumn.get_num_aligned()
                    sna[curr_position] = num_align
                    # print(x.get_num_aligned())
                    ##### Average Base Quality #####
                    # TOOD: Graph consensus base qualities
                    #  all base qualities for position, take average
                    pos[curr_position] = curr_position
                    if num_align > 0:
                        bq[curr_position] = sum(pileColumn.get_mapping_qualities()) / num_align
                    ##### Errors at position #####
                    pos_indel[curr_position] = self.mindel(pileColumn)
                elif curr_ref != pileColumn.reference_name:
                    ##### store values for sequence pileup ########
                    self.sarr[curr_ref] = sna
                    self.baq[curr_ref] = bq
                    self.consensus_pos[curr_ref] = pos
                    self.indel[curr_ref] = pos_indel
                    self.curr_ref = pileColumn.reference_name
                    ##### clear tmp arrays ######
                    ref_len = self.bamfile.header.get_reference_length(curr_ref)
                    sna = [0] * ref_len
                    bq = [0] * ref_len
                    pos = [-1] * ref_len
                    pos_indel = [-1] * ref_len
                    i=0
            except StopIteration:
                ##### store values for sequence pileup ########
                self.sarr[curr_ref] = sna
                self.baq[curr_ref] = bq
                self.consensus_pos[curr_ref] = pos
                self.indel[curr_ref] = pos_indel
                self.curr_ref = pileColumn.reference_name
                ##### clear tmp arrays ######
                ref_len = self.bamfile.header.get_reference_length(curr_ref)
                sna = [0] * ref_len
                bq = [0] * ref_len
                pos = [-1] * ref_len
                pos_indel = [-1] * ref_len
                i = 0
                break

        data["sarr"] = self.sarr
        data["baq"] = self.baq
        data["consensus_pos"] = self.consensus_pos
        data["indel"] = self.indel
        return data

    #### Count Modifications at ALL Bases at Query Position and Notation for isMatch/isInsert #####
    def mindel(self, x):
        mindel = x.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
        #print("Mindel")
        #print(mindel)
        deletion = 0
        insertion = 0
        mismatch = 0
        mismatches = []
        for b in mindel:
            # print(b)
            if re.search('\*|-', b) != None:
                deletion += 1
                # print("Deletion")
            elif re.search('\+', b) != None:
                insertion += 1
                # print("insertion")
            elif re.search('[ACGTacgt]', b) != None:
                mismatch += 1
                mismatches.append(b)
        ####### TODO: will need all mismatches for correlation
        return deletion, insertion, mismatch, Counter(mismatches).most_common(1)

    ########### Print Summary of Modifications for each Reference Sequence #################
    def print_summary(self):
        samfile = pysam.AlignmentFile(self.f_bam, 'rb')
        for key, item in self.sarr.items():
            # reads per sequence
            reads = samfile.fetch(key)
            print("Reads: " + str(len(list(reads))))
            # print reference sequence length
            print("Sequence: " + key)
            print("Reference Sequence Length: " + str(self.bamfile.header.get_reference_length(key)))
            # print reference position of modification
            print("Consensus Alignment Positions (-1 no reads aligned to position): " + str(self.consensus_pos[key]))
            # print modifications as tuple
            print("Modifications by position: Deletion/Insertion/Mismatch (-1 no reads aligned to position): " \
                  + str(self.indel[key]))
            print("Total number of alignments by position: " + str(item))


    ######### Modification Rate per base for all reads ###############
    ## Output: Average Modification rate per sequence by nucleotide,
    ##         Summary plots of modfication by position, and average modification
    def print_summary_modification_rates(self, plot=False):
        base_header = []
        ave_mod_header = ['Condition','Insertion', 'Deletion', 'Mismatch', 'Average']

        #print("Printing Modification Rates...")
        mod_del= {}
        mod_mis = {}
        mod_ins = {}
        mod_base = {}
        ave_mod_base = {}
        for seq_name, mods in self.indel.items():
            base_header.append(seq_name)
            seq_len  = self.bamfile.header.get_reference_length(seq_name)
            #print(key + ": " + str(seq_len))
            # temp base mod rates
            overall_mod = []
            ins_mod = []
            del_mod = []
            mis_mod = []
            for mod_tuple, num_aligned in zip(mods, self.sarr[seq_name]):
                if num_aligned > 0:
                    overall_mod.append((sum(mod_tuple[0:3]) / num_aligned))
                    ins_mod.append((mod_tuple[0] / num_aligned))
                    del_mod.append((mod_tuple[1] / num_aligned))
                    mis_mod.append((mod_tuple[2] / num_aligned))
                else:
                    ##### num_aligned == 0, indicated with -1 no reads, no errors indicated with 0 #####
                    overall_mod.append(-1)
                    ins_mod.append(-1)
                    del_mod.append(-1)
                    mis_mod.append(-1)
                    #print(0)
            ###### store modification rates by sequence and position #########
            mod_base[seq_name] = overall_mod
            mod_del[seq_name] = del_mod
            mod_mis[seq_name] = mis_mod
            mod_ins[seq_name] = ins_mod
            ######### Average Modification Rates by type of modification #################
            i1 = statistics.mean([ elem for i,elem in enumerate(ins_mod) if elem != -1])
            d1 = statistics.mean([ elem for i,elem in enumerate(del_mod) if elem != -1])
            m1 = statistics.mean([ elem for i,elem in enumerate(mis_mod) if elem != -1])
            o1 = statistics.mean([i1,d1,m1])
            ave_mod_base[seq_name] = [self.dir_name, i1, d1,m1,o1]

        ###### Create Output Folder ifne ##########
        save_path = self.save_path + self.dir_name + '_Modification_Rates/'
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        ##### Print CSV ###########
        #TODO: remove outliers (1) from averages in plots
        df = pd.DataFrame.from_dict(ave_mod_base, orient="index", columns=ave_mod_header)
        df.to_csv(save_path + self.dir_name + "_ave_mod_rates.csv")
        df2 = pd.DataFrame()
        for key, values in mod_base.items():
            dft = pd.DataFrame()
            dft['position'] = np.arange(0,len(values))
            dft['contig'] = key
            dft['Basecall_Reactivity'] = values
            dft['Quality'] = self.baq[key]
            areads = np.array(self.sarr[key])
            x = np.array(mod_mis[key])
            dft['Mismatch'] = x #np.where((areads > 0), x / areads, x)
            x = np.array(mod_del[key])
            dft['Deletion'] = x #np.where(areads > 0, x / areads, x)
            x = np.array(mod_ins[key])
            dft['Insertion'] = x #np.where(areads > 0, x / areads, x)
            dft['Aligned_Reads'] = areads
            dft['Sequence'] = [char for char in self.reference_sequences[key]]
            df2 = pd.concat([df2, dft])

        df2.to_csv(save_path + self.dir_name + "_mod_rates.csv")
        df2 = pd.DataFrame.from_dict(mod_base, orient="index")
        df2.unstack()
        tf = df2.T

        if plot is True:
            pltmod.plot_average_mod_rate(df, dir_name=self.dir_name, save_path=self.save_path)
            pltmod.plot_average_mod_by_pos_rate(tf, dir_name=self.dir_name, save_path=self.save_path)

    ## Input: list of consensus positions, list of modifications by nucleotides
    ## indel[sequence_name] = (deletions, insertions, mismatches), modification, save path
    ## Output: Plots by sequence and nucleotide of each type of modification (optional)
    ## Note: positions with no reads are excluded
    def print_all_modification_rates(self,plot=False):
        position = self.consensus_pos
        indel = self.indel

        deletions = []
        insertions = []
        mismatch = []
        mismatches = []
        summary = []

        ###### Create Output Folder ifne ##########
        save_path = self.save_path + self.dir_name + '_Modification_Rates/'
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        for key in position:  # select current sequence
            sequence = [char for char in self.reference_sequences[key]]  # get nucleotides
            ### for each position get errors ###
            for i in range(len(position[key])):
                if (type(indel[key][i]) is int) or (indel[key][i]==-1):
                    dl = -1
                    ins = -1
                    mis = -1
                    mnum = -1
                    called_base = -1
                else:
                    dl = indel[key][i][0]
                    ins = indel[key][i][1]
                    mis = indel[key][i][2]
                    mnum = (0 if not indel[key][i][3] else indel[key][i][3][0][1])
                    called_base = (sequence[i] if not indel[key][i][3] else re.sub("\W", "", indel[key][i][3][0][0]))
                deletions.append((position[key][i],dl))
                insertions.append((position[key][i],ins))
                mismatch.append((position[key][i],mis))
                # pos, number of mismatches,mismatch nucleotide,refseq nucleotide
                mismatches.append(
                    (position[key][i],mnum,called_base,sequence[i]))

                summary.append((position[key][i], dl, ins, mis))
                # if indel[key][i] != -1:  ##### base not aligned
                #     if indel[key][i][0] > 0:
                #         deletions.append((position[key][i], indel[key][i][0]))
                #     if indel[key][i][1] > 0:
                #         insertions.append((position[key][i], indel[key][i][1]))
                #     if indel[key][i][2] > 0:
                #         mismatch.append((position[key][i], indel[key][i][2]))
                #     if len(indel[key][i][3]) > 0:
                #         # pos, number of mismatches,mismatch nucleotide,refseq nucleotide
                #         mismatches.append(
                #             (position[key][i], indel[key][i][3][0][1], re.sub("\W", "", indel[key][i][3][0][0]), sequence[i]))
                #     if (indel[key][i][0] > 0 or indel[key][i][1] > 0 or indel[key][i][2] > 0):
                #         summary.append((position[key][i], indel[key][i][0], indel[key][i][1], indel[key][i][2]))
            ### if modifications then plot ########
            if len(deletions) > 0:
                df = pd.DataFrame(deletions, columns=["Position", "Deletions"])
                df.to_csv(save_path + key + "_deletion_rates.csv")
                if plot is True:
                    pltmod.plot_modification(deletions, seq_name=key, mod_type="Deletions", dir_name=self.dir_name,
                                         save_path=self.save_path)
            if len(insertions) > 0:
                df = pd.DataFrame(insertions, columns=["Position", "Insertion Rate"])
                df.to_csv(save_path + key + "_insertion_rates.csv")
                if plot is True:
                    pltmod.plot_modification(insertions, seq_name=key, mod_type="Insertions",
                                         dir_name=self.dir_name, save_path=self.save_path)
            if len(mismatch) > 0:
                df = pd.DataFrame(mismatch, columns=["Position", "Mismatch Rate"])
                df.to_csv(save_path + key + "_mismatch_rates.csv")
                if plot is True:
                    pltmod.plot_modification(mismatch, seq_name=key, mod_type="Mismatch",
                                         dir_name=self.dir_name, save_path=self.save_path)
            if len(mismatches) > 0:
                #### calculate mismatch bias ####
                # mismatch_bias(mismatches) ### In use
                df = pd.DataFrame(mismatches, columns=["Position", "Number of Mismatches", "Modified Nt", "Reference Nt"])
                df.to_csv(save_path + key + "_mismatch_bias_rates.csv")
                if plot is True:
                    pltmod.plot_mismatch(mismatches, seq_name=key, dir_name=self.dir_name,
                                     save_path=self.save_path)
            ##### Plot Summary of Modifications ############
            if len(summary) > 0:
                df = pd.DataFrame(summary, columns=["Position", "Deletion Rate", "Insertion Rate", "Mismatch Rate" ])
                df.to_csv(save_path + key + "_summary_rates.csv")
                if plot is True:
                    pltmod.plot_modification_summary(summary, dir_name=self.dir_name, seq_name=key,
                                                 save_path=self.save_path)
            deletions = []
            insertions = []
            mismatch = []
            mismatches = []
            summary = []
