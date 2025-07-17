import os
import re
import mysql.connector
import mariadb
import sys
import DashML.Database_fx.DB as db
import pandas as pd
import datetime
import numpy as np
import uuid
from sqlalchemy import text
import Insert_DB as dbins




dtr = pd.DataFrame.from_dict(
                {"contig":  ['test'],
                "sequence": ['ACGTU'],
                "sequence_name": ['test_sequence'],
                "sequence_len": 5,
                "temp": 37,
                "is_modified": 0,
                "type1": ['dms09'],
                "type2": ['dms010'],
                "complex": 0,
                "run": 1}, orient='columns')

def insert_orphan_structure():
        con = db.get_pd_remote()
        ts = datetime.datetime.now()
        sequence ='GGATCTTCGGGGCAGGGTGAAATTCCCGACCGGTGGTATAGTCCACGAAAGTATTTGCTTTGATTTGGTGAAATTCCAAAACCGACAGTAGAGTCTGGATGAGAGAAGATTC'
        sequence_len= len(sequence)
        with (con.connect() as conn):
            with conn.begin():
                    #insert df of sid, psotion, contig, sequence, and experiment=putative
                    nt = list(sequence.strip().replace("'", ""))
                    position = np.arange(0, sequence_len)
                    df_structure = pd.DataFrame({'SID': 12,
                                                 'contig': 'FMN',
                                                 'position': position, 'sequence': nt, 'experiment': 'putative',
                                                 'timestamp':ts})
                    #print(df_structure)
                    df_structure.to_sql('Structure', if_exists='append', index=False, con=conn)
                    conn.commit()


insert_orphan_structure()
# print(dtr)
# dtr['ID'] = 84
# print(dtr['ID'][0].astype(str))
# contig_value = (dtr["ID"].astype(str) + " " + dtr["contig"].astype(str) + " "
#                             + dtr["type1"].astype(str) + " " + dtr["type2"].astype(str))[0]
# print(contig_value)
# sys.exit(0)
#res = dbins.insert_library(dtr)
#print(res)
