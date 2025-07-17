import re
import sys
import DashML.Database_fx.DB as db
import pandas as pd
from sqlalchemy import text


def select_unmod_ssi():
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_UnmodSSI ()'
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_library_full():
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = "CALL Select_LibraryFull();"
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        raise Exception(error)
    finally:
        conn.close()
        return df


def select_library(contig="HCV", temp=37,type1="dmso", type2="acim", complex=0):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectLibrary ("{0}", {1}, "{2}", "{3}", {4})'.format(contig,temp,type1, type2, complex)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df

def select_librarybyid(lid):
    try:
        con = db.get_pd_remote()
        lid = re.sub(r"[\[\]\s\n'\"']", '', str(lid))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectLibraryById("{0}")'.format(lid)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df


def select_lids(temp=37,type1="dmso", type2="acim", complex=0, contig= "HCV"):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectLibraryLids ("{0}", {1}, "{2}", "{3}", {4})'.format(contig,temp,type1, type2, complex)
            df = pd.read_sql_query(query, conn)
            conn.commit()


    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_predict(lids):
    try:
        con = db.get_pd_remote()
        lids = str(lids).replace("[", "").replace("]", "").replace(' ', '')

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectPredict ("{0}")'.format(lids)
            print(query)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        return None
    finally:
        conn.close()
        return df


def select_unmod(seq, temp=37, type1='dmso', type2='dmso', complex=0):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectUnmod ("{0}",{1},"{2}","{3}", "{4}")'.format(seq,temp,type1,type2,complex)
            #print(query)
            df = pd.read_sql_query(query, conn)
            print(seq + " records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_unmod_lid(lids):
    try:
        con = db.get_pd_remote()
        lids = str(lids).replace("[", "").replace("]", "").replace(' ', '')

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_UnmodLID ("{0}")'.format(lids)
            #print(query)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute select unmod by lid: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_mod_lid(lids):
    try:
        con = db.get_pd_remote()
        lids = str(lids).replace("[", "").replace("]", "").replace(' ', '')

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ModLID ("{0}")'.format(lids)
            #print(query)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute select mod by lid: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df

def select_mod(seq, temp=37, type1='dmso', type2='acim', complex=0):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectMod ("{0}",{1},"{2}","{3}", "{4}")'.format(seq,temp,type1,type2,complex)
            #print(query)
            df = pd.read_sql_query(query, conn)
            print(seq + " records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_bc_mod(lids):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectBC ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)
            print("Basecall Mod records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


#### Grouped Unmodified Values ######
def select_bc_unmod(lids):
#grouped data
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectBCUnmod ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)
            print("Basecall Unmod records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df



def select_peaks(mod_lids, unmod_lids):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectPeaks ("{0}", "{1}")'.format(mod_lids, unmod_lids)
            df = pd.read_sql_query(query, conn)
            print("Peaks records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_peaks_unmod(mod_lids, unmod_lids):
    try:
        con = db.get_pd_remote()

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_PeaksUnmod ("{0}", "{1}")'.format(mod_lids, unmod_lids)
            df = pd.read_sql_query(query, conn)
            print("Peaks records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df

#Select Average Predictions from Read Depth
def select_read_depth_ave(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDAverage("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df

#Select Reads that are average length or longer
def select_read_depth_full(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectReadDepthFull ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)
            print("Read Depth Full records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_read_depth_full_all(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDFLID ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df

def select_max_clusters(lids):
    try:
        con = db.get_pd_remote()
        lids =  re.sub('\[|\]|\\s|\\n+|\'', '', lids)

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_MaxClusters ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df

##### Select Centroids with Control Structures #######
def select_centroids(lids):
    try:
        con = db.get_pd_remote()

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectCentroids ("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)
            print("Centroid records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


##### Select Centroids with Reactivity #######
def select_centroidz(lids):
    try:
        con = db.get_pd_remote()
        lids =  re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_Centroids ("{0}")'.format(lids)
            print(query)
            df = pd.read_sql_query(query, conn)
            print("Centroid records " + str(len(df)))

    except Exception as error:
        print("Failed to execute {0}: {1}".format(query, error))
    finally:
        conn.close()
        return df

def select_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', lids)

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectStructure("{0}")'.format(lids)
            print(query)
            df = pd.read_sql_query(query, conn)
            print("Structure records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_readdepth(lids):
    try:
        con = db.get_pd_remote()

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ReadDepth("{0}")'.format(lids)
            print(query)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute read depth: {}".format(error))
        return None
    finally:
        conn.close()
        return df


def select_secondarystructures(lids=None):
    try:
        con = db.get_pd_remote()

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            if lids is None: # select all
                query = 'CALL Select_SecondaryStructures()'
                df = pd.read_sql_query(query, conn)
                print("Structure Secondary records " + str(len(df)))
            else: #select list
                lids = re.sub('\[|\]|\\s|\\n+|\'\(\)', '', str(lids))
                query = 'CALL Select_SecondaryStructure({0})'.format(lids)
                df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df


def select_max_structure_probabilities(ssid):
    try:
        con = db.get_pd_remote()


        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectMaxStructureProbabilities({0})'.format(ssid)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df



def select_continued_reads(lids, lids2=None):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        if lids2 is None:
            lids2 = lids
        else:
            lids2 = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids2))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RdfContinue("{0}", "{1}")'.format(lids, lids2)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df


def select_putativestructures(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            query = 'CALL Select_PutativeStructures({0})'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute {0}: {1}".format(query,error))
    finally:
        conn.close()
        return df

#### VARNA output format
# dvarna['position'] = dvarna['position'] + 1
# dvarna.to_csv(save_path + seq + "_VARNA.csv", index=False, header=False)
def select_varna(seq_ids):
    try:
        con = db.get_pd_remote()
        seq_ids = str(seq_ids.insert(0, " "))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectPredict ("{0}")'.format(seq_ids)
            df = pd.read_sql_query(query, conn)
            print("Select Predict records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df



##### FROM TABLE READ_DEPTH ####
#### Vienna .dat output format
# dvienna = dfr[['position', 'RNAFold_Shape_Reactivity']]
# dvienna['position'] = dvienna['position'] + 1
# dvienna.to_csv(save_path + seq + "_rnafold.dat", index=False, sep='\t', header=False)
def select_vienna(seq_ids):
    try:
        con = db.get_pd_remote()
        seq_ids = str(seq_ids.insert(0, " "))

        with con.connect() as conn:
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL SelectPredict ("{0}")'.format(seq_ids)
            df = pd.read_sql_query(query, conn)
            print("Vienna records " + str(len(df)))

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
        return df
