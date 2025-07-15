import sys, re
import platform
import os.path
import traceback
import pandas as pd
from platformdirs import user_documents_path
from pathlib import Path

######### Parse Reference Sequence from FASTA ###############
def parse_fasta(filepath):
    sequence = {}
    sequence_name = ""
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                s = re.split('[\s>]', line)
                sequence_name = s[1]
                # print(s[1])
            elif re.search("(?i)^[ACGTU]", line):
                sequence[sequence_name] = line[:-1].upper()
                # print(line)
    return sequence

# reverse sequence direction for nanopolish signal analysis
def parse_fasta_signal(filepath):
    sequence = {}
    sequence_name = ""
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                s = re.split('[\s>]', line)
                sequence_name = s[1]
                # print(s[1])
            elif re.search("(?i)^[ACGTU]", line):
                sequence[sequence_name] = line[:-1].upper()[::-1]
    return sequence

######### Alignment Paths #####################

def check_file(f_path, f_name, ext):
    for filename in os.listdir(f_path):
        if filename is f_name:
            return f_name, True
        if filename.endswith(ext):
            return filename, True
    return None, False
def set_path_error(f_path, modification):
    ######### Path and Directory format ##############
    f_path = f_path
    print(user_documents_path())
    output = user_documents_path() / "DTLandscape_Output"
    output.mkdir(parents=True, exist_ok=True)
    subdir = output / "Alignment_Out"
    subdir.mkdir(parents=True, exist_ok=True)
    f_save_path = str(subdir.resolve()) + '/'
    #verify fasta
    f_name, _ = check_file(f_path, "library.fasta", '.fasta')
    if f_name == None:
        raise Exception("Fasta file not found!")
    else:
        f_reference = os.path.join(f_path, f_name)

    #verify sam
    f_name, _ = check_file(f_path, "aln.sam", '.sam')
    if f_name == None:
        raise Exception("Sam file not found!")
    else:
        f_sam = os.path.join(f_path, f_name)

    #verify baq
    f_name, _ = check_file(f_path, "aln.baq.bam", '.baq.bam')
    if f_name == None:
        raise Exception("BAQ.BAM file not found!")
    else:
        f_bam = os.path.join(f_path, f_name)

    reference_sequences = parse_fasta(f_reference)
    dir_name = modification.upper()
    return reference_sequences, f_path, f_reference, f_sam, f_bam, dir_name, f_save_path
