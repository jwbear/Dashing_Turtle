import sys
import os
import pandas as pd
import platform
from platformdirs import user_documents_path
from pathlib import Path
from DashML.Basecall.Basecalls import *
from DashML.Basecall import Basecall_Paths as dp
import DashML.Database_fx.Select_DB as dbsel
import DashML.Database_fx.Insert_DB as dbins


######### Set file paths and output names ##############
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

##### Get Base Modifications ######
# f_path: path to Alignment folder for sequence, containing alignments
# modification  : chemical modification type (any string)
# plot : whether to plot alignment plots
def getModifications(lid, contig, f_path, modification="DMSO", plot=False):
    reference_sequences, f_path, f_reference, f_sam, f_bam, dir_name, save_path = \
        dp.set_path_error(f_path=f_path, modification=modification)

    #### check contig in reference_sequences
    if reference_sequences.get(contig) is None:
        raise Exception(f"Sequence name '{contig}' not listed in library FASTA file. Please check library and try again")
        return None

    print(f_path, modification, "Running basecall....")

    plot = plot

    #================== Get Basecall Modifications/Errors =====================
    #### input: name of dataset for datat_files to retrieve, name of directory with specific modifications ####
    ### all reads in fasta are aligned by default ###
    pa = ParseAlignment(reference_sequences, f_path, f_reference, f_sam, f_bam, dir_name, save_path)
    #outputs single molecule analysis
    print("Getting modifications by read.......")
    pa.get_base_modifications_by_read()
    ## get base modifications for all reference sequences in fasta ##
    ## returns dict of modifications for each sequence ##
    print("Getting summary base modifications.......")
    data = pa.get_base_modifications()
    ## nice print summary of all base modfications for all sequences in fasta ##
    print("Printing summary.......")
    pa.print_summary()
    ## print summary of overall and average base modfication rates for all sequences in fasta to csv##
    ## optionally plot rates ##
    print("Plotting summary data.......")
    pa.print_summary_modification_rates(plot=plot)
    ## print all base modfication rates for all sequences in fasta to csv##
    ## optionally plot rates to png ##
    pa.print_all_modification_rates(plot=plot)
    #insert <Modification>_mod_rates into basecall data
    subdir = Path(save_path) / (modification + "_Modification_Rates")
    subdir.mkdir(parents=True, exist_ok=True)
    f = modification +"_mod_rates.csv"
    f = subdir / f"{f}"
    f= str(f.resolve())
    #print("fp", f)
    df = pd.read_csv(f)
    df = df.loc[df['contig'] == contig]
    if len(df) <= 0 :
        raise Exception(contig + ' not found in alignment file. Check that the sequence name/ contig matches the '
                                 'library (fasta/fa) sequence name. Try adding the sequence again.')
        return None
    else:
        df = df[['position','contig','Basecall_Reactivity','Quality','Mismatch',
                 'Deletion','Insertion','Aligned_Reads','Sequence']]
        df['LID'] = lid
        dbins.insert_basecall_rates(df)
        return save_path


#### single sequence processing ####
# Input: full path to directory containing library.fasta, aln.sam, aln.baq.bam
# multiple contigs in a fastq select only matching lid+contig pair

def check_lids(lids):
    return dbins.check_lids(lids)
def get_modification(lid, contig, f_path, modification, plot):
    try:
        save_path = None
        is_lid = check_lids(lid)
        if is_lid <= 0:
            raise Exception("Sequences must be in library before adding basecall data...")
        save_path = getModifications(lid=lid, contig=contig, f_path=f_path, modification= modification, plot=plot)
        if save_path is not None:
            plot1 = save_path + modification + "_Modification_Plots/" + "Average_Modification_Rate.png"
            plot2 = save_path + modification + "_Modification_Plots/" + "Position_Modification_Rate.png"
            return plot1, plot2
        else:
            raise Exception("Basecall Data Insufficient")
            return None
    except Exception as err:
        raise Exception('Basecall Error: ' + str(err))
        return None


#get_modification(52, contig="FMN", f_path='/Users/timshel/structure_landscapes/DashML/Basecall/DMSO/fmn/Alignment',
#              modification='DMSO', plot=True)
# sys.exit(0)
