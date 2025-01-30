import pandas as pd

# Defining a header to be used for all files.
HEADER = f"#AA\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10\ts11\tS12\n"
# Defining columns for pandas dataframe.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]
numCols = len(COL_NAMES)

# Simple function to write to file.
# Parameters:
#   df = The data frame to be written.
#   fileName = A string representing the name of the file to be written to.
def writeDf(df, fileName):
    outputFile = open(file = fileName, mode = "w")
    outputFile.write(HEADER)
    outputFile.write(df.to_csv(sep = '\t', header = False))
    outputFile.close()

# This function will take a dictionary and print out all of its information
# into separate files. I am aiming to make it output a computationally friendly
# file format.
# Parameters:
#   sectorDict = The dictionary containing sectorContacts objects.
#   dictType = The name of the dictionary (bb2bb, bb2sc, etc.).
def printSectorInfo(sectorDict, dictName):
    for parCode in sectorDict.keys():
        # Initializing a data frame with no rows.
        countDf = pd.DataFrame(columns = COL_NAMES)
        # Gathering information.
        sectorContacts_obj = sectorDict[parCode]
        sectorContacts_sectors = sectorContacts_obj.sectors
        # Building dataframe.
        for sector in sectorContacts_sectors.keys():
            sectorAAs = sectorContacts_sectors[sector]
            for AA in sectorAAs.keys():
                contactCount = sectorAAs[AA]
                try:
                    countDf.loc[AA, sector] = contactCount
                except:
                    countDf.loc[AA] = [0] * numCols
                    countDf.loc[AA, sector] = contactCount
        # Making a file names.
        countFileName = f"{parCode}.{dictName}.counts.txt"
        sectorProbFileName = f"{parCode}.{dictName}.sectorProbs.txt"
        AAProbFileName = f"{parCode}.{dictName}.AAProbs.txt"
        # Last preparations for writing data frames.
        countDf = countDf.fillna(0)
        sectorProbDf = countDf.div(countDf.sum(axis = 0), axis = 1) # Probability distribution for each sector
        sectorProbDf = sectorProbDf.fillna(0)
        AAProbDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for each contacting AA
        AAProbDf = AAProbDf.fillna(0)
        # Writing frames.
        writeDf(countDf, countFileName)
        writeDf(sectorProbDf, sectorProbFileName)
        writeDf(AAProbDf, AAProbFileName)