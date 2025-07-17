import sys, os, re
import subprocess
import pandas as pd
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel


# RNAeval -v -d2 < input.txt
# Calulate MFE of native psuedoknot free structures in Structure_Secondary Table

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


def get_MFE(df):
    try:
        lid = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['LID'].unique()))
        id = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['ID'].unique()))
        seq = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['contig'].unique()))
        print(seq)
        sequence = re.sub('\[|\]|\\s|\\n+|\'|\'', '', str(df['sequence'].unique()))
        structure = re.sub('\[|\]|\\s|\\n+|\'', '', str(df['secondary'].unique()))
        line = sequence + "\n" + structure + "\n"
        print(line)
        #send to rnaeval
        p1 = subprocess.Popen(["RNAeval", "-v", "-d2"],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        #reval = sequence + "\n" + ss + "\n@"
        p1.communicate(line.encode())[0]
        #p1.communicate(input=reval.encode())
        output = str(p1.communicate(timeout=300)).split("\\n")
        p1.stdout.close()
        p1.stdin.close()
        #save mfes to db
        mfe = extract_mfe(structure,output)
        ## insert mfe
        dbins.insert_secondary_mfe(id, lid, mfe)
    except subprocess.CalledProcessError as e:
        print(e)
        return None
    except Exception as e:
        print(e)
        raise Exception(str(e))
        return None
    finally:
        p1.stdout.close()
        p1.stdin.close()

def extract_mfe(ss,output):
    # get sequences and lengths
    for line in output:
        if ss in line:
            l = line.strip().split(' ')
            mfe = re.sub('\(|\)', '', l[1])
            print(mfe)
            break
    return mfe


#### get native mfes for secondary structures
def get_mfes(lids=None):
    if lids is None: #select all
        df = dbsel.select_secondarystructures()
        for id in df['ID'].unique():
            mfe_t = df['mfe'].loc[df['ID'] == id][0]
            if mfe_t is None:
                print(id)
                get_MFE(df[df['ID']==id])
    else:
        df = dbsel.select_secondarystructures(lids)
        for id in df['ID'].unique():
            mfe_t = df['mfe'].loc[df['ID'] == id][0]
            if mfe_t is None:
                print(id)
                get_MFE(df[df['ID']==id])

# get_mfes('51')
# sys.exit(0)
