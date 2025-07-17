import os
import re
import numpy as np
import mysql.connector
import mariadb
import sys
import DashML.Database_fx.DB as db
import pandas as pd
import datetime
import uuid
from sqlalchemy import text


def check_lids(lids):
    try:
        con = db.get_pd_remote()
        lids = str(lids).replace("[", "").replace("]", "")

        with con.connect() as conn:
            #check if exists in library
            with conn.begin():
                query = 'CALL SelectLibraryById("{0}")'.format(lids)
                table = conn.execute(text(query)).fetchall()
                if len(table) <= 0:
                    raise Exception("Sequences must be in library before adding basecall data...")
                conn.commit()

    except Exception as error:
        print("Failed to execute insert basecall error rates {}: {}".format(query, error))
    finally:
        conn.close()
        return len(table)

def insert_library(df):
    try:
        last_insert_id = None
        con = db.get_pd_remote()
        ts = datetime.datetime.now()
        df['timestamp'] = ts
        #`contig`, `sequence`, `sequence_name`, `sequence_len`, `temp`, `complex`, `is_modified`, `is_putative
        contig =  str(df['contig'].unique()).replace("[", "").replace("]", "")
        sequence =  str(df['sequence'].unique()).replace("[", "").replace("]", "")
        secondary = str(df['secondary'].unique()).replace("[", "").replace("]", "")
        experiment = str(df['experiment'].unique()).replace("[", "").replace("]", "")
        if (experiment == 'None') or (experiment == ''):
            experiment == 'putative'
        sequence_name =  str(df['sequence_name'].unique()).replace("[", "").replace("]", "")
        type1 = str(df['type1'].unique()).replace("[", "").replace("]", "")
        type2 = str(df['type2'].unique()).replace("[", "").replace("]", "")
        sequence_len =  int(str(df['sequence_len'].unique()).replace("[", "").replace("]", ""))
        temp = float(re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['temp'].unique())))
        run = int(re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['run'].unique())))
        complex =  0
        is_modified = int(re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['is_modified'].unique())))
        is_putative = 1


        with (con.connect() as conn):
            with conn.begin():
                #insert to library
                query = ('CALL Ins_Library ({0},{1},{2},{3},{4},{5},{6},{7},{8}, {9},"{10}", {11}, @lid)'). \
                         format(contig, sequence, sequence_name, sequence_len,temp,type1, type2, complex,
                                is_modified,is_putative,datetime.datetime.now(), run)
                last_insert_id = conn.execute(text(query)).first()
                last_insert_id = int(last_insert_id[0])
                #print("ID", last_insert_id)
                if last_insert_id==None:
                    raise Exception("No primary key returned.")
                else:
                    #insert secondary structure
                    if secondary not in ('None', ''):
                        df_secondary = pd.DataFrame.from_dict({'LID': [last_insert_id],
                                                     'contig': [sequence_name.strip().replace("'", "")],
                                                     'secondary': [secondary],
                                                     'timestamp': [ts]} , orient='columns')
                        df_secondary.to_sql('Structure_Secondary', if_exists='append', index=False, con=conn)
                    #insert to library joiner table
                    query = ('CALL Ins_StructureLibrary ({0},{1}, @msid)').format(last_insert_id, contig)
                    structureid = conn.execute(text(query)).first()
                    structureid = int(structureid[0])
                    #print("SID", structureid)
                    #print(query)
                    if structureid==None:
                        raise Exception("No structure id returned.")
                    else:
                        #insert df of sid, psotion, contig, sequence, and experiment=putative
                        nt = list(sequence.strip().replace("'", ""))
                        seqlen = len(nt)
                        position = np.arange(0, sequence_len)
                        df_structure = pd.DataFrame.from_dict({'SID': np.full(seqlen,structureid),
                                                     'contig': np.full(seqlen, contig.strip().replace("'", "")),
                                                     'position': position,
                                                     'sequence': nt,
                                                     'experiment': np.full(seqlen, experiment.strip().replace("'", "")),
                                                     'timestamp': np.full(seqlen,ts)},
                                                      orient = 'columns')
                        #print(df_structure)
                        df_structure.to_sql('Structure', if_exists='append', index=False, con=conn)

                        conn.commit()

    except Exception as error:
        print("Failed to execute insert library: {}".format(error))
    finally:
        return last_insert_id
        conn.close()

def insert_signal(df):
    try:
        #check columns
        cols = ['LID', 'contig', 'position', 'reference_kmer', 'read_index', 'event_level_mean', 'event_length',
                'event_stdv']
        dfcols = df.columns
        for c in cols:
            if c not in dfcols:
                raise Exception("Signal file missing col: " + c)

        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        #split df into unmod/mod
        df_unmod = df[df['type1']==df['type2']]
        df_unmod.drop(columns=['type1', 'type2'], inplace=True)
        df_mod = df[df['type1'] != df['type2']]
        df_mod.drop(columns=['type1', 'type2'], inplace=True)
        lids = re.sub(r"[\[\]\s\n'\"']", '', str(df['LID'].unique()))

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_Signal("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            with conn.begin():
                #insert signal, into mod/unmods based on lids, combine delete
                if len(df_unmod) > 0:
                    df_unmod.to_sql('Unmodified', if_exists='append', index=False, con=conn)
                if len(df_mod) > 0:
                    df_mod.to_sql('Modified', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute insert signal: {}".format(error))
    finally:
        conn.close()
        return

# Column names must match with alchemy
def insert_rx_full(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_ReactivityFull("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            with conn.begin():
                df.to_sql('Reactivity_full', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Reactivity_full WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_read_depth_full_clear(df, continue_reads=False):
    ## long running separate function
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")
        if 'base_pair_prob' not in df.columns:
            df['base_pair_prob'] = 0

        with con.connect() as conn:
            # remove old data
            if continue_reads == False:
                with conn.begin():
                    query = 'CALL Del_ReadDepthFull("{0}")'.format(lids)
                    conn.execute(text(query))
                    conn.commit()

    except Exception as error:
        raise Exception("Failed to execute insert read depth: {}".format(error))
    finally:
        conn.close()
        return


#update existing rdf Predict after running BPP

def insert_read_depth_full_update(lids, threshold, rx_threshold=.7):
    try:
        con = db.get_pd_remote()
        lids = str(lids).replace("[", "").replace("]", "")

        with con.connect() as conn:
            with conn.begin():
                query = ('CALL UPDATE_BPP({0},{1},{2});').format(lids, threshold, rx_threshold)
                print(query)
                conn.execute(text(query))
                conn.commit()

    except Exception as error:
        print("Failed to execute insert read depth probability update: {}".format(error))
    finally:
        conn.close()

def insert_read_depth_full(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        #lids = str(df['LID'].unique()).replace("[", "").replace("]", "")
        if 'base_pair_prob' not in df.columns:
            df['base_pair_prob'] = 0

        with con.connect() as conn:
            with conn.begin():
                df.to_sql('Read_depth_full', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute insert read depth: {}".format(error))
        return
    finally:
        conn.close()
        return


def insert_read_depth(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]","")

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_ReadDepth("{0}")'.format(lids)
                conn.execute(text(query))
                df.to_sql('Read_depth', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_basecall_rates(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = re.sub(r"[\[\]\s\n'\"']", '', str(df['LID'].unique()))

        with con.connect() as conn:
            #check if exists in library
            with conn.begin():
                query = 'CALL SelectLibraryById("{0}")'.format(lids)
                table = conn.execute(text(query)).fetchall()
                if len(table) <= 0:
                    raise Exception("Sequences must be in library before adding basecall data...")
                conn.commit()

            #remove old data
            with conn.begin():
                query = 'CALL Del_Basecall("{0}")'.format(lids)
                conn.execute(text(query))
                df.to_sql('Basecall', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute insert basecall error rates {}: {}".format(query, error))
    finally:
        return

#### seqids format "[,1,2,3,4]"
def insert_basecall(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")

        with con.connect() as conn:
            #remove old data
            with conn.begin():
                query = 'CALL Del_BasecallPeaks("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            #insert new data
            with conn.begin():
                df.to_sql('Basecall_Peaks', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Basecall_Peaks WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()



def insert_peaks(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")

        with con.connect() as conn:
            #remove old data
            with conn.begin():
                query = 'CALL Del_Peaks("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            #insert new data
            with conn.begin():
                df.to_sql('Peaks', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry ="SELECT * FROM Peaks WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_lof(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")


        with con.connect() as conn:
            #remove old data
            with conn.begin():
                query = 'CALL Del_Lof("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            #insert new data
            with conn.begin():
                df.to_sql('Lof', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Lof WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()




def insert_gmm(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = str(df['LID'].unique()).replace("[", "").replace("]", "")


        with con.connect() as conn:
            #remove old data
            with conn.begin():
                query = 'CALL Del_Gmm("{0}")'.format(lids)
                conn.execute(text(query))
                conn.commit()

            #insert new data
            with conn.begin():
                df.to_sql('Gmm', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Gmm WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_clusters(df, seq_ids):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        df = df.drop(columns=['read_id'])
        method = str(df['method'].unique()).replace("[", "").replace("]","").replace(' ','')
        seq_ids = str(seq_ids).replace("[", "").replace("]", "").replace(' ','')

        with con.connect() as conn:
            #remove old data
            with conn.begin():
                query = 'CALL Del_Clusters("{0}", {1})'.format(seq_ids, method)
                conn.execute(text(query))
                conn.commit()

            #insert new data
            with conn.begin():
                df.to_sql('Clusters', if_exists='append', index=False, con=conn)
                conn.commit()

            del df
            # display all records
            qry = 'SELECT * FROM Clusters WHERE FIND_IN_SET(LID,' + seq_ids + ') > 0;'
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_centroids(df, seq_ids):
    try:
        if len(df) <= 0:
            return
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        method = str(df['method'].unique()).replace("[", "").replace("]","").replace(' ','')
        seq_ids = str(seq_ids).replace("[", "").replace("]","").replace(' ','')

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_Centroids("{0}", {1})'.format(seq_ids, method)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Centroids', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()



def insert_centroid_distance(df):
    try:
        if len(df) <= 0:
            return

        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        clust_method = str(df['method'].unique()).replace("[", "").replace("]", "").replace(' ','')
        lids = list(df['LID'].unique())
        lids = str(lids).replace("[", "").replace("]", "").replace(' ','')
        clusters = list(df['cluster'].unique())
        clusters = str(clusters).replace("[", "").replace("]", "").replace(' ','')

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_CentroidDistance("{0}", "{1}", {2})'.format(lids, clusters, clust_method)
                #print(query)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Centroid_Distance', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Centroid_Distance  WHERE FIND_IN_SET(LID,'" + lids + "') > 0;"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

            del df

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()



def insert_secondary_mfe(id, lid, mfe):
    try:
        con = db.get_pd_remote()
        timestamp = datetime.datetime.now()

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Ins_Mfe({0}, {1}, {2}, "{3}")'.format(id, lid, mfe, timestamp)
                conn.execute(text(query))
                conn.commit()

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()


def insert_centroid_secondary(df):
    try:
        if len(df)<=0:
            return

        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        cluster_method = str(df['method'].unique()).replace("[", "").replace("]", "").replace(' ', '')
        lids = list(df['LID'].unique())
        lids = str(lids).replace("[", "(").replace("]", ")").replace(' ', '')
        clusters = list(df['cluster'].unique())
        clusters = str(clusters).replace("[", "(").replace("]", ")").replace(' ', '')

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_CentroidSecondary({0}, {1}, {2})'.format(lids, clusters, cluster_method)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Centroid_Secondary', if_exists='append', index=False, con=conn)
                conn.commit()

            # display all records
            qry = "SELECT * FROM Centroid_Secondary  WHERE LID IN (" + lids + ");"
            table = conn.execute(text(qry)).fetchall()
            print(len(table))

            del df

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()


def insert_structure_secondary_intrx(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lid = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['LID'].unique()))
        lid2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['LID2'].unique()))
        read = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['read_index'].unique()))

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_StructureSecondaryIntrx({0}, {1}, {2})'.format(lid, lid2, read)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Structure_Secondary_Interactions', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute structure secondary: {}".format(error))

    finally:
        conn.close()


def insert_centroid_secondary_intrx(df):
    try:
        if len(df)<=0:
            raise Exception('Database empty', len(df))
            return

        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lid = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['LID1'].unique()))
        lid2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['LID2'].unique()))
        cluster = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['cluster1'].unique()))
        cluster2 = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['cluster2'].unique()))
        method = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['method'].unique()))

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = 'CALL Del_CentroidSecondaryIntrx({0}, {1}, {2}, {3})'.format(lid, lid2, cluster,
                                                                                     cluster2, method)
                #print(query)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Centroid_Secondary_Interactions', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        raise Exception("Failed to execute centroid secondary intrx: {}".format(error))

    finally:
        conn.close()



def insert_structure_bpp(df_bpp, df_prob):
    try:
        con = db.get_pd_remote()
        df_bpp['timestamp'] = datetime.datetime.now()
        df_prob['timestamp'] = datetime.datetime.now()
        # ssid = "0x" + re.sub('(UUID)|\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df_bpp['SSID'].unique()))
        # an_integer = int(ssid, 16)
        # ssid = hex(an_integer)

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                #Delete from structure_secondary_interactions to prevent orphans
                #query = "CALL Del_StructureBPP({0})".format(ssid)
                #print(query)
                #conn.execute(text(query))
                df_bpp.to_sql('Structure_BPP', if_exists='append', index=False, con=conn)
                conn.commit()

            with conn.begin():
                df_prob.to_sql('Structure_Probabilities', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute structure_bpp: {}".format(error))
    finally:
        conn.close()


def insert_centroid_bpp(df_bpp, df_prob):
    try:
        con = db.get_pd_remote()
        df_bpp['timestamp'] = datetime.datetime.now()
        df_prob['timestamp'] = datetime.datetime.now()
        # ssid = "0x" + re.sub('(UUID)|\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df_bpp['SSID'].unique()))
        # an_integer = int(ssid, 16)
        # ssid = hex(an_integer)

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                #Delete from structure_secondary_interactions to prevent orphans
                #query = "CALL Del_StructureBPP({0})".format(ssid)
                #print(query)
                #conn.execute(text(query))
                df_bpp.to_sql('Centroid_BPP', if_exists='append', index=False, con=conn)
                conn.commit()

            with conn.begin():
                df_prob.to_sql('Centroid_Probabilities', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute structure_bpp: {}".format(error))
    finally:
        conn.close()


def insert_centroid_regions(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()
        lids = re.sub('\[|\]|\\s|\\n+|\|\"|\'', '', str(df['LID'].unique()))
        method = re.sub('\(|\)|\[|\]|\\s|\\n+|\|\"|\'', '', str(df['method'].unique()))

        with con.connect() as conn:
            # remove old data
            with conn.begin():
                query = "CALL Del_CentroidRegions({0}, '{1}')".format(lids, method)
                #print(query)
                conn.execute(text(query))
                conn.commit()

            # insert new data
            with conn.begin():
                df.to_sql('Centroid_Regions', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute centroid regions: {}".format(error))

    finally:
        conn.close()


def insert_metric(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()

        with con.connect() as conn:
            with conn.begin():
                df.to_sql('Structure_Metric', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()



def insert_secondary_metric(df):
    try:
        con = db.get_pd_remote()
        df['timestamp'] = datetime.datetime.now()

        with con.connect() as conn:
            with conn.begin():
                df.to_sql('Structure_Metric_Secondary', if_exists='append', index=False, con=conn)
                conn.commit()

    except Exception as error:
        print("Failed to execute: {}".format(error))
        sys.exit(0)
    finally:
        conn.close()
