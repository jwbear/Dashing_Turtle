import sys, os, re, ast
import subprocess
import pandas as pd
import numpy as np
import uuid
from snowflake import SnowflakeGenerator
import time
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel


##### Calculate Base Pairing Probabilities for Predictions either single molecule reactivities (read index = -1)
# or no reactivity molecules (DMSO read index = -2)  #####
####### USED BY Predict.py


#RNA basepairing with reactivity is probably better here
# RNAfold -p -d2 --noLP --MEA --shape=HCV_rnafold2.dat < hcv.fa > hcv_bp.out
# RNAcofold -a -d2 --noLP < sequences.fa > cofold.out
# bp percentages where predict is true but over 95% can be unmodified, may be due to cluster effects
# ignore non-predicted or missing values
#RNAfold -p -d2 --noLP < test_sequenc.fa > test_sequenc.out
# RNAcofold -a -d2 --noLP < sequences.fa > cofold.out
# $ RNAfold --shape=reactivities.dat < sequence.fa
# where the file reactivities.dat is a two column text file with sequence positions (1-based)
# normalized reactivity values (usually between 0 and 2. Missing values may be left out, or assigned a negative score:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

gen = SnowflakeGenerator(42)

def get_tmp_dir():
    package_dir = os.path.dirname(os.path.abspath(__file__))  # Folder containing this file
    tmp_dir = os.path.join(package_dir, "TMP")
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir

#extract putative mfes from output
def extract_mfes(output, dt, dt2):
    contig = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(dt['contig'].unique()))
    lid = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(dt['ID'].unique()))
    ssid = re.sub('(UUID)|\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(dt['SSID'].unique()))
    read = int(re.sub('\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(dt['read_index'].unique())))

    contig2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(dt2['contig'].unique()))
    lid2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(dt2['ID'].unique()))
    ssid2 = re.sub('(UUID)|\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(dt2['SSID'].unique()))

    df = pd.DataFrame(columns=['SSID', 'LID', 'contig', 'LID2', 'contig2', 'read_index', 'secondary', 'mfe',
                               'frequency', 'deltag', 'type', 'mea'])
    # Define the conversion dictionary
    convert_dict = {'LID': int, 'contig': str, 'LID2': int, 'contig2': str, 'read_index': int,
                    'secondary': str, 'mfe' : float,'frequency' : float, 'deltag' : float, 'type' : str, 'mea': float}

    # Convert columns using the dictionary
    df = df.astype(convert_dict)
    mfe, freq, deltag, mea = 0, 0, 0, 0
    ab, aa, bb, a, b = "", "","","",""
    for i, line in enumerate(output.split("\\n")):
        #print(line)
        secondary, type  = "", ""
        if re.search('(frequency)', line):
            l = re.split(";", line,1)
            #print(l)
            f = re.search('([0-9]+.[0-9]+e*-*[0-9]*)', l[0])
            freq = f.group()
            #print("frequency", freq)
            d = re.search('-*([0-9]+\.[0-9]+e*-*[0-9]*)', l[1])
            deltag = d.group()
            #print("deltag", deltag)
        elif re.search(r'(-*[0-9]+\.[0-9]+\\t)', line):
            l =  re.split(r'\\t', line)
            ab, aa, bb, a, b = l[0], l[1], l[2], l[3], l[4]
            #print("energies", ab, aa, bb, a, b)
        elif re.search('^((\.+)|(\(+|\)+))\W', line):
            l = re.split("\s", line,1)
            secondary = l[0]
            #print("structure", secondary)
            energy = l[1]
            if re.search('^(\(|-[0-9])', energy):
                type = "MFE"
                mfe = re.sub('\(|\)', '', energy)
                #print("mfe", mfe)
            elif re.search('^(\[|-[0-9])',  energy):
                type = "PROB"
                mfe = re.sub('\[|\]','', energy)
                #print("prob mfe", mfe)
            elif re.search('^(\{|-[0-9])',energy):
                i = 0
                e = re.split("\s", energy)
                if len(e[0]) <= 1:
                    i = 1
                if re.search('(d\=)', energy):
                    type = "CENTROID"
                    mfe = re.sub('\{|\}|d|\=', '', e[i])
                    #print("centroid mfe",mfe)
                    distance = re.sub('\{|\}|d|\=', '', e[i+1])
                    #print("centroid distance", distance)
                elif re.search('(MEA\=)', energy):
                    type = "MEA"
                    mfe = re.sub('\{|\}|(MEA)|\=', '', e[i])
                    #print("mea mfe", mfe)
                    mea = re.sub('\{|\}|(MEA)|\=', '', e[i + 1])
                    #print("MEA", mea)
        #'LID', 'contig', 'secondary', 'mfe', 'cluster', 'frequency', 'diversity', 'type', 'distance'
        if secondary != "":
            if len(mfe) == 0:
                mfe = 0
            df.loc[len(df)] = [ssid, lid, contig, lid2, contig2, read, secondary, mfe, freq, deltag, type, mea]

    df.loc[df['type'] == 'MFE', 'frequency'] = float(freq)
    df.loc[df['type'] == 'MFE', 'deltag'] = float(deltag)
    #print(df)

    #insert into db
    dbins.insert_structure_secondary_intrx(df)
    return {'AB': ab, 'AA': aa, 'BB' : bb, 'A': a, 'B': b}

def get_bpfiles(contig, contig2, read_index):
    file_list = {
        'AA' : "AA{0},{1}{2}_dp5.ps".format(contig, contig2, read_index),
        'AB' :"AB{0},{1}{2}_dp5.ps".format(contig, contig2, read_index),
        'BB' : "BB{0},{1}{2}_dp5.ps".format(contig, contig2, read_index),
        'A' : "A{0},{1}{2}_dp5.ps".format(contig, contig2, read_index),
        'B' : "B{0},{1}{2}_dp5.ps".format(contig, contig2, read_index)
    }
    return file_list

#extract base pairing probabilities from output save to Structure_BPP
#  ubox The upper right triangle displays all predicted base pairs with not more than two inconsistent sequences.
# lbox The lower left triangle contains only the secondary structure formed by the most believable base pairs.
# It still contains pairs that are not in the final prediction, lbox all equal
def extract_bpps(df, df2, mfes):
    contig = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['contig'].unique()))
    ssid = re.sub('(UUID)|\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df['SSID'].unique()))
    seqlen = int(re.sub('\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['sequence_len'].unique())))
    read_index = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['read_index'].unique()))

    contig2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['contig'].unique()))
    seqlen2 = int(re.sub('\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['sequence_len'].unique())))


    # get file list by creating file names
    bp_files = get_bpfiles(contig, contig2, read_index)

    # for each file loop through and extract bpp and probabilities
    for rx, bpf in bp_files.items():
        #print(rx)
        mfe = mfes[rx]
        base1, base2, prob = [], [], []
        def get_bases():
            try:
                tmp_dir = get_tmp_dir()
                tmp_file_path = os.path.join(tmp_dir, bpf)

                with open(tmp_file_path, "r+") as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if re.search("([0-9]+ [0-9]+ [0-9]+\.*[0-9]* (ubox))", line.strip()):
                        for l in lines[i:]:
                            if re.search("(showpage)|(end)|(%%EOF)", l):
                                return
                            else:
                                bases = l.strip().split()
                                base1.append(int(bases[0]))
                                base2.append(int(bases[1]))
                                prob.append(float(bases[2]))
                                #print(bases)
            except Exception as e:
                print(e)
                raise Exception(str(e))
                return
            finally:
                f.close()

        get_bases()
        bppid = next(gen)
        df_bpp = pd.DataFrame({'BPPID': bppid, 'SSID': ssid, 'type_interaction': rx, 'mfe':float(mfe)}, index=[0])
        df_prob = pd.DataFrame({'BPPID': bppid, 'SSID': ssid, 'base1':base1, 'base2':base2, 'probability':prob})
        # fix bases RNAcofold extends base numbers for different molecules, AA
        df_prob['base1'] = np.where(df_prob['base1'] > seqlen, df_prob['base1'] - seqlen, df_prob['base1'])
        df_prob['base2'] = np.where(df_prob['base2'] > seqlen2, df_prob['base2'] - seqlen2, df_prob['base2'])
        #insert into db
        dbins.insert_structure_bpp(df_bpp, df_prob)
    return

##### get intra and inter molecular interactions
def get_dimers(df, df2):
    try:
        #create input file
        ## vienna rna does not always follow naming conventions if names exceed length
        tmp = get_tmp_dir()

        contig = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['contig'].unique()))
        sequence = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['sequence'].unique()))
        read_index = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['read_index'].unique()))


        contig2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['contig'].unique()))
        sequence2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['sequence'].unique()))

        if len(sequence + sequence2) >= 5000:
            print('Length Exceeds Vienna CoFold')
            # TODO switch to RNAFold only to calculate putative structure, no base pairing info
            # cenfl limit is 5k
            #RNAfold -p -d2 --noLP < sequence1.fa > sequence1.out
            return None

        # create input file
        sf = "{}{}{}{}_sequence.fa".format(tmp, contig, contig2, read_index)
        f = open(sf, "w")
        lines = "> {},{}{}\n{}&{}\n".format(contig, contig2, read_index, sequence, sequence2)
        f.write(lines)
        f.close()

        #send to rnacofold
        in_path = sf
        #print(in_path)
        #out_path = save_path_dmso + seqname + "_" + seqname + "_intra.out"
        p1 = subprocess.Popen(["RNAcofold", "-p", "-a", "-d2", "--noLP", "--MEA", in_path],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=tmp)
        output = str(p1.communicate(timeout=None))
        #print(output)
        # generate unique run number
        ssid = next(gen)
        df['SSID'] = ssid
        df2['SSID'] = ssid
        # extract mfes and base pairing data
        mfes = extract_mfes(output, df, df2)
        # mfes = {'AB': .5, 'AA': .7, 'BB': .8, 'A': .9, 'B': .11}
        extract_bpps(mfes=mfes, df=df, df2=df2)
    except subprocess.CalledProcessError as e:
        print(e)
        return None
    except Exception as e:
        print(e)
        sys.exit(1)
    finally:
        if 'f' in locals():
            f.close()
        if 'p1' in locals():
            p1.stdout.close()
            p1.stdin.close()


def get_reactive_dimers(df, df2, reactivity):
    try:
        # create input file
        ## vienna rna does not always follow naming conventions if names exceed length
        tmp = get_tmp_dir()

        contig = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['contig'].unique()))
        sequence = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['sequence'].unique()))
        read_index = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['read_index'].unique()))

        contig2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['contig'].unique()))
        sequence2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df2['sequence'].unique()))

        if len(sequence + sequence2) >= 5000:
            print('Length Exceeds Vienna CoFold')
            #TODO pesky switch to RNAFold only to calculate putative structure, no base pairing info
            # limit is 5k nt
            #RNAfold -p -d2 --noLP < sequence1.fa > sequence1.out
            return None

        # create input file
        tmp_dir = get_tmp_dir()
        filename = f"{contig}{contig2}{read_index}_sequence.fa"
        sf = os.path.join(tmp_dir, filename)
        f = open(sf, "w")
        lines = "> {},{}{}\n{}&{}\n".format(contig, contig2, read_index, sequence, sequence2)
        f.write(lines)
        f.close()

        # create dat file of reactivities
        tmp_dir = get_tmp_dir()
        filename = f"{contig}{contig2}{read_index}.dat"
        datf = os.path.join(tmp_dir, filename)
        f = open(datf, "w")
        for i in range(0, len(reactivity)):
            f.write(str(i + 1) + "\t" + str(reactivity[i]) + "\n")
        f.close()
        dat_file = "--shape=" + os.path.abspath(datf)
        print("dat file: ", dat_file)

        # send to rnacofold
        in_path = sf
        # print(in_path)
        # out_path = save_path_dmso + seqname + "_" + seqname + "_intra.out"
        p1 = subprocess.Popen(["RNAcofold", "-p", "-a", "-d2", "--noLP", "--MEA", dat_file, in_path],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=tmp)
        output = str(p1.communicate(timeout=None))
        #print(output)
        # generate unique run number
        ssid = next(gen)
        df['SSID'] = ssid
        df2['SSID'] = ssid
        # extract mfes and base pairing data
        mfes = extract_mfes(output, df, df2)
        # mfes = {'AB': .5, 'AA': .7, 'BB': .8, 'A': .9, 'B': .11}
        extract_bpps(mfes=mfes, df=df, df2=df2)
        return ssid

    except subprocess.CalledProcessError as e:
        print(e)
        return None
    except Exception as e:
        print(e)
        return None
    finally:
        if 'f' in locals():
            f.close()
        if 'p1' in locals():
            p1.stdout.close()
            p1.stdin.close()


def get_probabilities(lid, lid2, reactivity=None, read=-2):
    # select molecules
    df = dbsel.select_librarybyid(lid)
    df['read_index'] = read
    df2 = dbsel.select_librarybyid(lid2)
    df2['read_index'] = read
    #suffix = str(lid) + str(lid2) + str(np.abs(read))
    if (len(df)<=0) or (len(df2) <=0):
        raise Exception("Sequence not found in library: LID {}, LID {}".format(lid,lid2))
        return None

    if reactivity is None:
        ssid = get_dimers(df, df2)
    else:
        ssid = get_reactive_dimers(df, df2, reactivity)
    return ssid

# sys.exit(0)
