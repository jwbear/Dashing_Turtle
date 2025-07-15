# import required modules
import os
import mysql.connector
import mariadb
import sys



# create connection object
# Connect to MariaDB Platform
def create_connection():
    try:
        # local
        con = mysql.connector.connect(
            user="root",
            password="",
            host="127.0.0.1",
            port=3306,
            database="DASH"
        )

        #remote
        # con = mysql.connector.connect(
        #     user="jwbear",
        #     password= "fU8VYQSv9!GGH9vM",
        #     host="127.0.0.1",
        #     port=28901,
        #     database="jwbear_DASH"
        # )

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
        path = '/Users/timshel/structure_landscapes/DashML/Deconvolution/BC/'
        for f in os.listdir(path):
            fname = path + f
            if ('metric' not in f) and ('bc_' not in f) and (f.endswith(".csv")):
                print(fname)
                querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Basecall " +
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
                querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Signal_ " +
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
                querybc = ("LOAD DATA LOCAL INFILE \"" + fname + "\" INTO TABLE DASH.Dwell " +
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
                querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Gmm " +
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
                querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Lof_dwell " +
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
                querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Lof_signal " +
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
        query2 = "select * from lof_dwell"
        cursor.execute(query2)

        # display all records
        table = cursor.fetchall()
        print(len(table))

        query2 = "select * from lof_signal"
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
        querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Structure " +
        "FIELDS TERMINATED BY ',' " +
        "LINES TERMINATED BY '\n' " +
        'IGNORE 1 ROWS ' +
        '(@ix,@position,@sequence,@basetype,@structuretype,@contig) '  +
        'SET ' +
        'position=@position,contig=@contig,sequence=@sequence,base_type=@basetype,' +
                                                'structure_type=@structuretype; \n')

        # executing cursor
        # NOTE: execute multi is not stable
        cursor.execute(querybc)
        con.commit()

        # assign data query
        query2 = "select * from Structure"
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
        querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Structure " +
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
        query2 = "select * from Structure"
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

        querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Shape " +
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

        querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Shape " +
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
        query2 = "select * from SHAPE"
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

import data_fx as dfx
# dfx.get_structure()
# dfx.get_structure_ext()
#dfx.get_shapemap()
#dfx.get_shape()
#dfx.get_shape_continuous()
#dfx.structuresforputativeseqs()

# load_putative_structure()
# load_shape()
# load_structure()
# load_lof()
# load_gmm()
# load_dwell()
# load_signal()
load_basecall()
sys.exit(0)
