import pandas as pd

# Defining a header to be used for all files.
HEADER = f"#AA\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10\ts11\tS12\n"
# Defining columns for pandas dataframe.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]
numCols = len(COL_NAMES)

# This function will take a dictionary and print out all of its information
# into separate files. I am aiming to make it output a computationally friendly
# file format.
# Parameters:
#   sectorDict = The dictionary containing sectorContacts objects.
#   dictType = The name of the dictionary (bb2bb, bb2sc, etc.).
def printSectorInfo(sectorDict, dictName):
    for parCode in sectorDict.keys():
        # Making a file name.
        fileName = f"{parCode}_{dictName}.txt"
        # Initializing a data frame with no rows.
        printDf = pd.DataFrame(columns = COL_NAMES)
        # Gathering information.
        sectorContacts_obj = sectorDict[parCode]
        sectorContacts_sectors = sectorContacts_obj.sectors
        # Building dataframe.
        for sector in sectorContacts_sectors.keys():
            sectorAAs = sectorContacts_sectors[sector]
            for AA in sectorAAs.keys():
                contactCount = sectorAAs[AA]
                try:
                    printDf.loc[AA, sector] = contactCount
                except:
                    printDf.loc[AA] = [0] * numCols
                    printDf.loc[AA, sector] = contactCount
        # Outputting to file.
        printDf = printDf.fillna(0)
        outputFile = open(file = fileName, mode = "w")
        outputFile.write(HEADER)
        outputFile.write(printDf.to_csv(sep = '\t', header = False))
        outputFile.close()