import sys, os, re
import platform
import subprocess
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

MAX_THREADS = 100

####
# Calculates putative structure and MFE of cluster centroids using reactivities
# Used to measure distance in MFE from native structure

#RNA basepairing with reactivity is probably better here
# RNAfold -p -d2 --noLP --MEA --shape=HCV_rnafold2.dat < hcv.fa > hcv_bp.out
# RNAcofold -a -d2 --noLP < sequences.fa > cofold.out
# todo bp percentages where predict is true but over 95% can be unmodified
# ignore non-predicted or missing values
#RNAfold -p -d2 --noLP < test_sequenc.fa > test_sequenc.out
# RNAcofold -a -d2 --noLP < sequences.fa > cofold.out
# $ RNAfold --shape=reactivities.dat < sequence.fa
# where the file reactivities.dat is a two column text file with sequence positions (1-based)
# normalized reactivity values (usually between 0 and 2. Missing values may be left out, or assigned a negative score:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def get_tmp_dir():
    package_dir = os.path.dirname(os.path.abspath(__file__))  # Folder containing this file
    tmp_dir = os.path.join(package_dir, "TMP")
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir

def scale_reactivities(reactivities):
    min = reactivities.min()
    max = reactivities.max()
    smin = 0
    smax = 2

    reactivities = ((reactivities - min) / (max - min)) * (smax - smin) + smin
    return reactivities

#extract putative mfes from output
def extract_mfes(lid, seq, cluster, clust_type,output):
    df = pd.DataFrame(columns=['LID', 'contig', 'secondary', 'mfe', 'cluster',
                               'frequency', 'diversity', 'type', 'distance', 'mea', 'method'])
    # Define the conversion dictionary
    convert_dict = {'LID': int, 'contig': str, 'secondary': str, 'mfe' : float, 'cluster' : int,
                    'frequency' : float, 'diversity' : float, 'type' : str, 'distance': float,
                    'mea': float, 'method': str}
    # Convert columns using the dictionary
    df = df.astype(convert_dict)

    for i, line in enumerate(output.split("\\n")):
        secondary, type  = "", ""
        mfe, distance, freq, diversity, distance, mea = 0, 0 ,0, 0, 0, 0
        if re.search('(frequency)', line):
            l = re.split(";", line,1)
            #print(l)
            f = re.search('([0-9]+.[0-9]+e*-*[0-9]*)', l[0])
            freq = f.group()
            #print("frequency", freq)
            d = re.search('([0-9]+.[0-9]+e*-*[0-9]*)', l[1])
            diversity = d.group()
            #print("diversity", diversity)
            df.loc[df['type']=='MFE', 'frequency'] = float(freq)
            df.loc[df['type'] == 'MFE', 'diversity'] = float(diversity)
        elif re.search('(\.+)|(\(+|\)+)\W', line):
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
            df.loc[len(df)] = [lid, seq, secondary, mfe, cluster, freq, diversity, type, distance, mea, clust_type]

    #insert into db
    dbins.insert_centroid_secondary(df)
    return

##### get putative structures for cluster centroids #####
# RNAfold -p -a -d2 --noLP --MEA  < sequence.fa > sequence_cofold.out
def getRNAfold(lid, seqname, sequence, clust_num, clust_rx, clust_type):
    try:

        # Get the package root directory dynamically
        package_dir = os.path.dirname(os.path.abspath(__file__))

        # Create (or ensure) TMP inside package
        tmp_dir = get_tmp_dir()

        # Compose full file path
        tmp_file_path = os.path.join(tmp_dir, "sequence.fa")

        #create input file
        f = open(tmp_file_path, "w")
        f.write(">" + seqname + "-" + str(clust_num) + "\n" + sequence + "\n")
        f.close()

        ### dat file, print reactivities to properly formatted dat file
        tmp_file_path = os.path.join(tmp_dir, "tmp.dat")
        f = open(tmp_file_path, "w")
        seqlen = len(sequence)
        for i in range(0, seqlen):
            f.write(str(i + 1) + "\t" + str(clust_rx[i]) + "\n")
        f.close()
        dat_file = "--shape=" + tmp_file_path
        #print("dat file: ", dat_file)

        #send to rnacofold
        #in_path = save_path_dir + "/sequence.fa"
        #print(in_path)
        #out_path = save_path_dir + '/' + seqname + "_" + str(clust_num) + ".out"
        p1 = subprocess.Popen(["RNAfold", "-p","-S 1.2", "-d2", "--noLP", "--MEA", dat_file, "sequence.fa"],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=get_tmp_dir())
        #reval = sequence + "\n" + ss + "\n@"
        #p1.communicate(input=reval.encode())
        output = str(p1.communicate(timeout=300))
        extract_mfes(lid, seqname, clust_num, clust_type,output)
    except subprocess.CalledProcessError as e:
        print(e)
        return None
    except Exception as e:
        print(e)
        raise Exception(str(e))
        return None
    finally:
        #f.close()
        p1.stdout.close()
        p1.stdin.close()


# for each of k reactivities
# get rnafold data with varying reactivities from clusters
def get_putative_structure(lids):
    # dataframes of centroids
    df = dbsel.select_centroidz(lids)
    clusters = df['cluster'].unique()

    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        for cluster in clusters:
            dt = df.loc[(df['cluster']==cluster) & (df['method']=='kmeans')]
            for li in dt['LID'].unique():
                seq_name = re.sub('\[|\]|\\s|\\n+|\'|\"', '', str(df['contig'].unique()))
                sequence = re.sub('\[|\]|\\s|\\n+|\'|\"', '', str(df['sequence'].unique()))
                reactivities = scale_reactivities(dt['centroid'].to_numpy())
                executor.submit(getRNAfold, li, seq_name, sequence, cluster, reactivities, 'kmeans')

        for cluster in clusters:
            dt = df.loc[(df['cluster']==cluster) & (df['method']=='hamming')]
            for li in dt['LID'].unique():
                seq_name = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['contig'].unique()))
                sequence = re.sub('\[|\]|\\s|\\n+|\'|\"', '', str(df['sequence'].unique()))
                dt.loc[dt['centroid'] == 1, 'centorid'] = 0
                dt.loc[dt['centroid'] == -1, 'centroid'] = 2
                reactivities = dt['centroid'].to_numpy()
                executor.submit(getRNAfold,li, seq_name, sequence, cluster, reactivities, 'hamming')



#### TODO check reactivities (kmeans) from centroid files
#get_putative_structure(37)
#sys.exit(0)
