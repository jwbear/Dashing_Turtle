import re
import DashML.Database_fx.DB as db
import pandas as pd


#### Get Average Metrics ######
def select_shape_control(lids, type, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        type = re.sub('\[|\]|\\s|\\n+|\'', '', str(type))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ShapeControl("{0}","{1}", "{2}")'.format(lids, type, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")


#### Get GMM Metrics
def select_gmm_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Gmm_Control("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_gmm_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Gmm_Structure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_bc_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_BcControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_bc_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_BcStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_signal_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_SignalControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_signal_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_SignalStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")


def select_dwell_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_DwellControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_dwell_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_DwellStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofs_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofsControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofs_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofsStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofd_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofdControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofd_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofdStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rdf_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDFControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rdf_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDFStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rd_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ReadDepthControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rd_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ReadDepthStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_vienna_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ViennaControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")
def select_vienna_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ViennaStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

#### Get Metrics By Read ######



def select_dwell_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_DwellControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_dwell_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_DwellStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofs_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofsControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofs_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofsStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofd_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofdControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_lofd_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_LofdStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rdf_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDFControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rdf_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_RDFStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rd_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ReadDepthControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_rd_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ReadDepthStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")

def select_vienna_control(lids, types):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))
        types = re.sub('\[|\]|\\s|\\n+|\'', '', str(types))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ViennaControl("{0}", "{1}")'.format(lids, types)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")
def select_vienna_structure(lids):
    try:
        con = db.get_pd_remote()
        lids = re.sub('\[|\]|\\s|\\n+|\'', '', str(lids))

        with (con.connect() as conn):
            # call stored proc with sequence as param
            # return result set convert to df
            query = 'CALL Select_ViennaStructure("{0}")'.format(lids)
            df = pd.read_sql_query(query, conn)

    except Exception as error:
        print("Failed to execute: {}".format(error))
    finally:
        conn.close()
        return df
        print("MySQL connection is closed")
