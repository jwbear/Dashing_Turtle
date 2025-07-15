# import required modules
import os
import mysql.connector
import mariadb
import sys
import datetime


# create connection object
# Connect to MariaDB Platform
def create_connection():
    try:
        # local
        # con = mysql.connector.connect(
        #     user="root",
        #     password="",
        #     host="127.0.0.1",
        #     port=3306,
        #     database="DASH"
        # )

        #remote
        con = mysql.connector.connect(
            user="jwbear",
            password= "fU8VYQSv9!GGH9vM",
            host="127.0.0.1",
            port=28901,
            database="jwbear_DASH"
        )

        return con

    except mysql.connector.Error as e:
        print(f"Error connecting to DB Platform: {e}")
        sys.exit(1)

def load_basecall():
    con = create_connection()
    cursor = con.cursor()
    #querybc = ""
    try:
        #LOAD BASECALL DATA
        #fname =  "/Users/timshel/structure_landscapes/DashML/Deconvolution/BC/cen_3'utr_complex_weightcompare.csv"
        path = '/Users/timshel/structure_landscapes/DashML/Basecall/'
        for f in os.listdir(path):
            fname = path + f
            complex = 0
            if 'complex' in f:
                complex = 1
            #position, contig, Basecall_Reactivity, Quality, Mismatch, Deletion, Insertion, Aligned_Reads, Sequence
            if f.endswith(".csv"):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Basecall " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@position,@contig,@reactivity,@quality,@mis,@del,@ins,@aligned,@sequence) '  +
                'SET ' +
                'position=@position,contig=@contig,basecall_reactivity=@reactivity,quality=@quality,mismatch=@mis,'
                                                              'deletion=@del,insertion=@ins,aligned_reads=@aligned,'
                                                              'sequence=@sequence,complex='+str(complex)+ '; \n')

                if complex==1:
                    print(querybc)
                    sys.exit(0)

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from Basecall"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")


def load_basecall_peaks():
    con = create_connection()
    cursor = con.cursor()
    #querybc = ""
    try:
        #LOAD BASECALL DATA
        #fname =  "/Users/timshel/structure_landscapes/DashML/Deconvolution/BC/cen_3'utr_complex_weightcompare.csv"
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/BC/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and ('bc_' not in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Basecall_Peaks " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@indexs,@position,@contig,@is_peak,@peak_height,@ins,@mis,@Predict,@varna) '  +
                'SET ' +
                'position=@position,contig=@contig,is_peak=@is_peak,peak_height=@peak_height,insertion=@ins,' +
                 'mismatch=@mis,predict=@Predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from Basecall_Peaks"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_signal():
    con = create_connection()
    cursor = con.cursor()
    #querybc = ""
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/Signal/'
        for f in os.listdir(path):
            fname = path + f
            print(fname)
            if ('metric' not in f) and (f.endswith(".csv")):
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Signal_ " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@read_index,@position,@contig,@predict,@varna) '  +
                'SET ' +
                'read_index=@read_index,position=@position,contig=@contig,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from Signal_"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_dwell():
    con = create_connection()
    cursor = con.cursor()
    #querybc = ""
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/Dwell/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Dwell " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@read_index,@position,@contig,@predict,@varna) '  +
                'SET ' +
                'read_index=@read_index,position=@position,contig=@contig,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from Dwell"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")


def load_gmm():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/GMM/'
        for f in os.listdir(path):
            fname = path + f
            print(fname)
            if f.endswith(".csv"):
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Gmm " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@ix,@position,@contig,@read_depth,@percent_modified,@predict) '  +
                'SET ' +
                'position=@position,contig=@contig,read_depth=@read_depth,percent_modified=@percent_modified,predict=@predict; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from Gmm"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_lof():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/LOF/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and ('dwell' in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Lof_dwell " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@ix,@position,@contig,@read_index,@predict,@varna) '  +
                'SET ' +
                'position=@position,contig=@contig,read_index=@read_index,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()
            elif ('metric' not in f) and ('signal' in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Lof_signal " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@ix,@position,@contig,@read_index,@predict,@varna) '  +
                'SET ' +
                'position=@position,contig=@contig,read_index=@read_index,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from jwbear_Dash.lof_dwell"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

        query2 = "select * from jwbear_Dash.lof_signal"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_lof_signal():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/LOF/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and ('signal' in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Lof_signal " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@ix,@position,@contig,@read_index,@predict,@varna) '  +
                'SET ' +
                'position=@position,contig=@contig,read_index=@read_index,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        query2 = "select * from jwbear_Dash.Lof_signal"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_lof_dwell():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/LOF/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and ('dwell' in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Lof_dwell " +
                "FIELDS TERMINATED BY ',' " +
                "LINES TERMINATED BY '\n' " +
                'IGNORE 1 ROWS ' +
                '(@ix,@position,@contig,@read_index,@predict,@varna) '  +
                'SET ' +
                'position=@position,contig=@contig,read_index=@read_index,predict=@predict,varna=@varna; \n')

                # executing cursor
                # NOTE: execute multi is not stable
                cursor.execute(querybc)
                con.commit()

        # assign data query
        query2 = "select * from jwbear_Dash.Lof_dwell"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_structure():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        fname = '/Users/timshel/structure_landscapes/DashML/Deconvolution/Structure/structures_ext.csv'

        print(fname)
        querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Structure " +
        "FIELDS TERMINATED BY ',' " +
        "LINES TERMINATED BY '\n' " +
        'IGNORE 1 ROWS ' +
        '(@ix,@position,@sequence,@basetype,@structuretype,@contig) '  +
        'SET ' +
        'position=@position,contig=@contig,sequence=@sequence,base_type=@basetype,' +
                                                'structure_type=@structuretype,annotation_type="retrospective"; \n')

        # executing cursor
        # NOTE: execute multi is not stable
        cursor.execute(querybc)
        con.commit()

        # assign data query
        query2 = "select * from jwbear_DASH.Structure"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_putative_structure():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        fname = '/Users/timshel/structure_landscapes/DashML/Deconvolution/Structure/putative_structures.csv'

        print(fname)
        querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Structure " +
        "FIELDS TERMINATED BY ',' " +
        "LINES TERMINATED BY '\n' " +
        'IGNORE 1 ROWS ' +
        '(@ix,@contig,@sequence,@position) '  +
        'SET ' +
        'position=@position,contig=@contig,sequence=@sequence,type="putative"; \n')

        # executing cursor
        # NOTE: execute multi is not stable
        cursor.execute(querybc)
        con.commit()

        # assign data query
        query2 = "select * from jwbear_DASH.Structure"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_shape():
    con = create_connection()
    cursor = con.cursor()
    try:
        #LOAD DATA
        fname = '/Users/timshel/structure_landscapes/DashML/Deconvolution/ShapeMap/shape_maps.csv'
        shape_type='\"map\"'
        print(fname)

        querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Shape " +
        "FIELDS TERMINATED BY ',' " +
        "LINES TERMINATED BY '\n' " +
        'IGNORE 1 ROWS ' +
        '(@ix,@position,@sequence,@reactivity,@contig,@predict) '  +
        'SET ' +
        'position=@position,contig=@contig,sequence=@sequence,reactivity=@reactivity,' +
                                                'predict=@predict,shape_type=' + shape_type + '; \n')

        # executing cursor
        # NOTE: execute multi is not stable
        cursor.execute(querybc)
        con.commit()

        fname = '/Users/timshel/structure_landscapes/DashML/Deconvolution/Shape/shapes.csv'
        shape_type = '\"ce\"'
        print(fname)

        querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE jwbear_DASH.Shape " +
                   "FIELDS TERMINATED BY ',' " +
                   "LINES TERMINATED BY '\n' " +
                   'IGNORE 1 ROWS ' +
                   '(@ix,@position,@sequence,@reactivity,@contig,@predict) ' +
                   'SET ' +
                   'position=@position,contig=@contig,sequence=@sequence,reactivity=@reactivity,' +
                   'predict=@predict,shape_type=' + shape_type + '; \n')

        # executing cursor
        # NOTE: execute multi is not stable
        cursor.execute(querybc)
        con.commit()

        # assign data query
        query2 = "select * from jwbear_Dash.Shape"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_acim(seqmap=None):
    con = create_connection()
    cursor = con.cursor()
    try:
        data_path ='/Users/timshel/structure_landscapes/DashML/Deconvolution/Decon/'
        if seqmap==None:
            seqmap = {
                "HCV": data_path + "ACIM/HCV_decon_signal_preprocess.csv",
                "cen_3'utr": data_path + "ACIM/cen_3'utr_decon_signal_preprocess.csv",
                "cen_3'utr_complex" : data_path + "ACIM/cen_3'utr_complex_decon_signal_preprocess.csv",
                "cen_FL": data_path + "ACIM/cen_FL_decon_signal_preprocess.csv",
                "cen_FL_complex" : data_path + "ACIM/cen_FL_complex_decon_signal_preprocess.csv",
                "ik2_3'utr": data_path + "ACIM/ik2_3'utr_decon_signal_preprocess.csv",
                "ik2_3'utr_complex": data_path + "ACIM/ik2_3'utr_complex_decon_signal_preprocess.csv",
                "ik2_FL": data_path + "ACIM/ik2_FL_decon_signal_preprocess.csv",
                "ik2_FL_complex" : data_path + "ACIM/ik2_FL_complex_decon_signal_preprocess.csv",
                "RNAse_P" : data_path + "ACIM/RNAse_P_decon_signal_preprocess.csv",
                "T_thermophila": data_path + "ACIM/T_thermophila_decon_signal_preprocess.csv",
                "HSP70_HSPA1A_37": data_path + "ACIM/HSP70_HSPA1A_37_decon_signal_preprocess.csv",
                "HSP70_HSPA1A_42": data_path + "ACIM/HSP70_HSPA1A_42A_decon_signal_preprocess.csv",
                "HSP70_HSPA1A_42_Run_2": data_path + "ACIM/HSP70_HSPA1A_42B_decon_signal_preprocess.csv"
            }

        # LOAD DATA
        for s,f in seqmap.items():
            print(s)
            print(f)
            _timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
            if ('42' in s):
                _temp = 42
            else:
                _temp = 37

            r = s.find('Run_')
            if r != -1:
                _run = s[r + 4:]
            else:
                _run = 1

            queryacim = ("LOAD DATA LOCAL INFILE \"" + f + "\" INTO TABLE jwbear_DASH.Acim " +
                       "FIELDS TERMINATED BY ',' " +
                       "LINES TERMINATED BY '\n' " +
                       'IGNORE 1 ROWS ' +
                       '(@ix,@contig,@position,@reference_kmer,@read_index,@event_level_mean,@event_length,@event_stdv) ' +
                       'SET ' +
                       'position=@position,contig=@contig,read_index=@read_index,reference_kmer=@reference_kmer,'
                       'event_level_mean=@event_level_mean,event_length=@event_length,'
                       'event_stdv=@event_stdv,temp='+str(_temp)+',run='+str(_run)+',timestamp="'+_timestamp+'"; \n')

            # executing cursor
            # NOTE: execute multi is not stable
            cursor.execute(queryacim)
            con.commit()

        # assign data query
        query2 = "select * from jwbear_DASH.Acim"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

def load_dmso(seqmap=None):
    con = create_connection()
    cursor = con.cursor()
    try:

        if seqmap == None:
            data_path ='/Users/timshel/structure_landscapes/DashML/Deconvolution/Decon/'
            seqmap = {
                "HCV": data_path + "DMSO/HCV_decon_signal_preprocess.csv",
                "cen_3'utr": data_path + "DMSO/cen_3'utr_decon_signal_preprocess.csv",
                "cen_3'utr_complex" : data_path + "DMSO/cen_3'utr_complex_decon_signal_preprocess.csv",
                "cen_FL": data_path + "DMSO/cen_FL_decon_signal_preprocess.csv",
                "cen_FL_complex" : data_path + "DMSO/cen_FL_complex_decon_signal_preprocess.csv",
                "ik2_3'utr": data_path + "DMSO/ik2_3'utr_decon_signal_preprocess.csv",
                "ik2_3'utr_complex": data_path + "DMSO/ik2_3'utr_complex_decon_signal_preprocess.csv",
                "ik2_FL": data_path + "DMSO/ik2_FL_decon_signal_preprocess.csv",
                "ik2_FL_complex" : data_path + "DMSO/ik2_FL_complex_decon_signal_preprocess.csv",
                "RNAse_P" : data_path + "DMSO/RNAse_P_decon_signal_preprocess.csv",
                "HSP70_HSPA1A_37": data_path + "DMSO/HSP70_HSPA1A_37_decon_signal_preprocess.csv",
                "HSP70_HSPA1A_42": data_path + "DMSO/HSP70_HSPA1A_42_decon_signal_preprocess.csv",
                "T_thermophila": data_path + "DMSO/T_thermophila_decon_signal_preprocess.csv"
            }

        # LOAD DATA
        for s,f in seqmap.items():
            print(s)
            print(f)
            _timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
            if ('42' in s):
                _temp = 42
            else:
                _temp = 37

            r = s.find('Run_')
            if r != -1:
                _run = s[r + 4:]
            else:
                _run = 1

            queryacim = ("LOAD DATA LOCAL INFILE \"" + f + "\" INTO TABLE jwbear_DASH.Dmso " +
                       "FIELDS TERMINATED BY ',' " +
                       "LINES TERMINATED BY '\n' " +
                       'IGNORE 1 ROWS ' +
                       '(@ix,@contig,@position,@reference_kmer,@read_index,@event_level_mean,@event_length,@event_stdv) ' +
                       'SET ' +
                       'position=@position,contig=@contig,read_index=@read_index,reference_kmer=@reference_kmer,'
                       'event_level_mean=@event_level_mean,event_length=@event_length,'
                       'event_stdv=@event_stdv,temp='+str(_temp)+',run='+str(_run)+',timestamp="'+_timestamp+'"; \n')

            # executing cursor
            # NOTE: execute multi is not stable
            cursor.execute(queryacim)
            con.commit()

        # assign data query
        query2 = "select * from jwbear_DASH.Dmso"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

    except mysql.connector.Error as error:
        print("Failed to execute: {}".format(error))
    finally:
        con.commit()
        cursor.close()
        con.close()
        print("MySQL connection is closed")

# import data_fx as dfx
# dfx.get_structure()
# dfx.get_structure_ext()
#dfx.get_shapemap()
#dfx.get_shape()
#dfx.get_shape_continuous()
#dfx.structuresforputativeseqs()

# load_putative_structure()
# load_shape()
# load_structure()
# load_putative_structure()
# load_lof()
# load_lof_dwell()
# load_lof_signal()
# load_gmm()
# load_dwell()
# load_signal()
load_basecall()
# load_acim()
# load_dmso()
sys.exit(0)
