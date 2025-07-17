### Generate Predictions #####

#### Basecall Input Data ####
# TODO: add to pipeline flow/ automate and make local copy
# compute canada file processing follow Dashing Turtle Documentation
# Get Compressed Fast5 Data  (submit_guppy.sh)
# Input: uncompressed fast5 files
# Output: compressed fast5 files
# 1. Run guppy on fast5 files to create compressed fast5 data for alignment

# Create Alignments (submit_aln2.sh)
# Input: library.fasta containing list of sequences, compressed fast5 files of sequences
# Output: .sam, .bam alignment of fast5 files with base quality data, single.fastq files with alignment data
# 2. Run minimap2 to create alignment data


### Basecall Data ####
# Input: .sam, .bam alignment files, library.fasta
# Output: average, position modification plots, spreadsheets of insertion, deletion, mismatch, summary,
#   and average rates for each sequence, list of all modification rates for all sequences, list of average
#   modificatioon rates for all sequences, modification for each base for all reads (native only),
#   file output format Alignment_Out/
#       <ChemicalModification>_Modification_<Rates|Plots|Reads>/
#       <ChemicalModification<_file_content>.csv
# 3. run_basecall.py to get list of modifications by read and summarized modifications by position


#### Basecall Reactivity #########
# Processed by Predictors/Peak_Analysis_BC
# Reactivity = Modified(Average Reactivity at Position) - Unmodified(Average reactivity of Position)
# Single Modification Reactivity = Modified(Average Reactivity at Position)(.5 reads) - Modified(Average reactivity of Position)(.5 reads)

### Basecalling Single Modification ####
# How do we define reactivity,
# initial thoughts just randomly split data into 2 parts and take the difference in reactivity
# splitting data at aln step, dividing single fastq
# Note: DMSO data should not be highly reactive
# Verify any correlations with structure in unmodified data
# TODO: Update if find better method
