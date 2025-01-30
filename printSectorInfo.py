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
#   dict = The dictionary containing dataframe objects.
#   dictType = The name of the dictionary (bb2bb, bb2sc, etc.).
def printSectorInfo(dict, dictName):
    for parKey in dict.keys():
        countDf = dict[parKey] # Dataframes stored in the dictionary have count values
        # Making a file names.
        countFileName = f"{parKey}.{dictName}.counts.txt"
        sectorProbFileName = f"{parKey}.{dictName}.sectorProbs.txt"
        AAProbFileName = f"{parKey}.{dictName}.AAProbs.txt"
        # Last preparations for writing data frames.
        sectorProbDf = countDf.div(countDf.sum(axis = 0), axis = 1) # Probability distribution for each sector
        AAProbDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for each contacting AA
        # Writing frames.
        writeDf(countDf, countFileName)
        writeDf(sectorProbDf, sectorProbFileName)
        writeDf(AAProbDf, AAProbFileName)