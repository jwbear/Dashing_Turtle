# import required modules
import os
import mysql.connector
import sys

# create connection object
# Connect to MariaDB Platform
def create_connection():
    try:
        con = mysql.connector.connect(
            user="root",
            password="",
            host="127.0.0.1",
            port=3306,
            database="DASH"
        )

        return con

    except mysql.connector.Error as e:
        print(f"Error connecting to DB Platform: {e}")
        sys.exit(1)



# # create cursor object
# cursor = con.cursor()
#
# # assign data query
# query1 = "desc test"
#
# # executing cursor
# cursor.execute(query1)
#
# # display all records
# table = cursor.fetchall()
#
# # describe table
# print('\n Table Description:')
# for attr in table:
#     print(attr)
#
# # assign data query
# query2 = "select * from test"
#
# # executing cursor
# cursor.execute(query2)
#
# # display all records
# table = cursor.fetchall()
#
# # fetch all columns
# print('\n Table Data:')
# for row in table:
#     print(row[0], end=" ")
#     print(row[1], end=" ")
#     print(row[2], end=" ")
#     print(row[3], end="\n")

# closing cursor connection
#cursor.close()

# closing connection object
#con.close()

# cursor = con.cursor()
# try:
#     #cursor.callproc('ts2')
#     # print results
#     # print("Stored Procedure")
#     # for result in cursor.stored_results():
#     #     print(result.fetchall())
#
#     #LOAD BASECALL DATA
#     fname =  "/Users/timshel/structure_landscapes/DashML/Deconvolution/BC/cen_3'utr_complex_weightcompare.csv"
#     querybc = ("LOAD DATA INFILE \"" + fname + "\" INTO TABLE DASH.Basecall_1 " +
#     "FIELDS TERMINATED BY ',' " +
#     "LINES TERMINATED BY '\n' " +
#     'IGNORE 1 ROWS ' +
#     '(@indexs,@position,@contig,@is_peak,@peak_height,@ins,@mis,@Predict,@varna) '  +
#     'SET ' +
#     'position=@position,contig=@contig,is_peak=@is_peak,peak_height=@peak_height,insertion=@ins,' +
#      'mismatch=@mis,predict=@Predict,varna=@varna;')
#     # executing cursor
#     cursor.execute(querybc)
#     con.commit()
#     # assign data query
#     query2 = "select * from Basecall_1"
#
#     # executing cursor
#     cursor.execute(query2)
#
#     # display all records
#     table = cursor.fetchall()
#
#     for r in table:
#         print(r)
#
#
#
#
#
#
# except mysql.connector.Error as error:
#     print("Failed to execute stored procedure: {}".format(error))
# finally:
#     cursor.close()
#     con.close()
#     print("MySQL connection is closed")
