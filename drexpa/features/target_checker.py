import sqlite3
import pandas as pd
def create_connection(db_file):
    """Create a database connection to the SQLite database given by the db_file
    :param db_file: database db_file
    :return: Connection object or None"""

    global con
    global cursor

    try:
        con = sqlite3.connect(db_file)
        cursor = con.cursor()
        print("Successfully connected to SQLite")
    except sqlite3.Error as e:
        print(e)
        print("ERROR")

def TargetsAndType(db_file, file_name, HGNC, type):
#Drugs that are small molecules and target a gene of interest.
    
    create_connection(db_file)
    drugList = []
    frames = []
    for row in cursor.execute("SELECT Drug.drugID FROM Drug NATURAL INNER JOIN Interaction WHERE HGNC = ? and type = ?", (HGNC, type)):
        s = str(row)
        drugList.append(s.strip("'(,)"))

    for i in drugList:
        cursor.execute("SELECT drugName, Drug.drugID, actionType, HGNC FROM Drug NATURAL INNER JOIN Interaction WHERE Drug.drugID = ?", (i, ))
        df = pd.DataFrame(cursor.fetchall(), columns=["drugName", "drugID", "actionType", "HGNC"])
        drugPanelTT(file_name, df)
    con.close()

def drugPanelTT(fileName, df):
    with open(fileName, 'a+') as f:
        f.seek(0)
        first_line = f.readline().rstrip('\n')
        if first_line != "#Name\tTarget":
            f.write("#Name\tTarget" + "\n") # Retrieve all targets for a given compound.
        targetList = []
        print("Target list: ", targetList)
        print(df)
        for ind in df.index:
            drugId = df["drugID"][ind]
            aType = df["actionType"][ind]
            if aType == "Inhibitor":
                    aT = "inhibits"
            targetList.append(df["HGNC"][ind])
        
        f.write(drugId + "\t" + aT + "\t" + "\t".join(sorted(set(targetList))) + "\n")
    f.close()    

def drug(db_file, file_name, drugName):
#Drug of interest
    
    create_connection(db_file)
    cursor.execute("SELECT * FROM Drug WHERE drugName = ? ", (drugName, ))
    df = pd.DataFrame(cursor.fetchall(), columns=["drugID", "drugName", "type"])
    drugPanel(file_name, df)
    con.close()

def drugPanel(fileName, df):
    with open(fileName, 'w') as f:
        f.write("#Name\tTarget" + "\n") # Retrieve all targets for a given compound.
        for ind in df.index:
            drugId = df["drugID"][ind]
            targetList = []
            cursor.execute("SELECT drugID, actionType, HGNC FROM Interaction WHERE drugID = ?", (drugId, ))
            df2 = pd.DataFrame(cursor.fetchall(), columns=["drugID", "actionType", "HGNC"])
            for i in df2.index:
                aType = df2["actionType"][i]
                if aType == "Inhibitor":
                    aT = "inhibits"
                targetList.append(df2["HGNC"][i])
            f.write(drugId + "\t" + aT + "\t" + "\t".join(targetList) + "\n")
            targetList.clear()
    f.close()

def targetProfile(db_file, file_name, id_list, kd_limit = None, ki_limit = None, ic50_limit = None):
    """Retrieve target profiles for a list of drugs with measured binding affinities below set limits and write them to a drug panel file. 
    :param db_file: database db_file, file_name: name for drugpanel file, id_list: list of drug IDs, kd_limit: Upper limit for kd value for binding affinity, ki_limit: Upper limit for ki value for binding affinity, ic50_limit: Upper limit for ic50 value for binding affinity. 
    """

    create_connection(db_file)
    frames = []
    for i in id_list:
        print(f'Processing drug ID: {i}')        
        if kd_limit != None:
            cursor.execute("SELECT DISTINCT Interaction.drugID, actionType, Interaction.HGNC FROM ((Interaction INNER JOIN MeasuredFor ON Interaction.DrugID = MeasuredFor.DrugID) INNER JOIN BindingAffinity ON MeasuredFor.BindingReactionID = BindingAffinity.BindingReactionID) WHERE Interaction.DrugID = ? and  Kd_max < ?", (i, kd_limit))
            df1 = pd.DataFrame(cursor.fetchall(), columns=["drugID", "actionType", "HGNC"])
            frames.append(df1)

        if ki_limit != None: 
            cursor.execute("SELECT DISTINCT Interaction.drugID, actionType, Interaction.HGNC FROM ((Interaction INNER JOIN MeasuredFor ON Interaction.DrugID = MeasuredFor.DrugID) INNER JOIN BindingAffinity ON MeasuredFor.BindingReactionID = BindingAffinity.BindingReactionID) WHERE Interaction.DrugID = ? and  Ki_max < ?", (i, ki_limit))
            df2 = pd.DataFrame(cursor.fetchall(), columns=["drugID", "actionType", "HGNC"])
            frames.append(df2)

        if ic50_limit != None: 
            cursor.execute("SELECT DISTINCT Interaction.drugID, actionType, Interaction.HGNC FROM ((Interaction INNER JOIN MeasuredFor ON Interaction.DrugID = MeasuredFor.DrugID) INNER JOIN BindingAffinity ON MeasuredFor.BindingReactionID = BindingAffinity.BindingReactionID) WHERE Interaction.DrugID = ? and  Kd_max < ?", (i, ic50_limit))
            df3 = pd.DataFrame(cursor.fetchall(), columns=["drugID", "actionType", "HGNC"])
            frames.append(df3)
        print(f'Targets for drug ID {i} processed: {frames}')

    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame(columns=["drugID", "actionType", "HGNC"]) # Empty dataframe with expected columns
        
    drugPanelTP(file_name, df)
    frames.clear()
    con.close()
    #return df

def drugPanelTP(fileName, df):
    with open(fileName, 'a+') as f:
        f.seek(0)
        first_line = f.readline().rstrip('\n')
        if first_line != "#Name\tTarget":
            f.write("#Name\tTarget" + "\n") # Retrieve all targets for a given compound.
        targetList = []
        try:
            for ind in df.index:
                drugId = df["drugID"][ind]
                aType = df["actionType"][ind]
                if aType == "Inhibitor":
                    aT = "inhibits"
                targetList.append(df["HGNC"][ind])
            f.write(drugId + "\t" + aT + "\t" + "\t".join(set(targetList)) + "\n")
        except:
            print("No interactions with binding affinity measured below the given limit value.")
    f.close()

def targetProfileBA(db_file, file_name, id_list, ic50_limit = None):
    """Retrieve target profiles for a list of drugs with measured binding affinities below set limits and write them to a drug panel file. 
    :param db_file: database db_file, file_name: name for drugpanel file, id_list: list of drug IDs, kd_limit: Upper limit for kd value for binding affinity, ki_limit: Upper limit for ki value for binding affinity, ic50_limit: Upper limit for ic50 value for binding affinity. 
    """
    create_connection(db_file)
    for i in id_list:
        # print(i)
        frames = []

        if ic50_limit: 
            cursor.execute("SELECT DISTINCT drugID, HGNC FROM MeasuredFor INNER JOIN BindingAffinity ON MeasuredFor.BindingReactionID = BindingAffinity.BindingReactionID WHERE DrugID = ? and  IC50_max <= ?", (i, ic50_limit))
            df3 = pd.DataFrame(cursor.fetchall(), columns=["drugID", "HGNC"])
            ic50_str = str(ic50_limit)
            df3["drugID"] = df3["drugID"] + "_" + ic50_str + "_" + "nM" # Adds the IC50 value at the end of the drug ID in the drug panel file
            frames.append(df3)
            
        df = pd.concat(frames, ignore_index=True)
        print(f'Number of targets for drug ID {i}: {len(df.index)}')
        # print(df)
        drugPanelBA(file_name, df)
        frames.clear()
        yield df

    con.close()

def drugPanelBA(fileName, df):
    with open(fileName, 'a+') as f:
        f.seek(0)
        first_line = f.readline().rstrip('\n')
        if first_line != "#Name\tAction\tTarget\tTarget2":
            f.write("#Name\tAction\tTarget\tTarget2" + "\n") 
        targetList = []
        try:
            for ind in df.index:
                drugId = df["drugID"][ind]
                aT = "inhibits" # Verify that the drug of interest is an inhibitor
                targetList.append(df["HGNC"][ind])
            f.write(drugId + "\t" + aT + "\t" + "\t".join(sorted(set(targetList))) + "\n")
        except:
            print("No interactions with binding affinity measured below the given limit value.")
        # print(targetList)
    f.close()

# TEST
# targetProfileBA('DrugTargetInteractionDB.db', 'test.txt', ['CHEMBL2041933', 'CHEMBL2103875'], ic50_limit = 10000)