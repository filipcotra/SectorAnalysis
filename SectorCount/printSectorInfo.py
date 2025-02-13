# Simple function to write to file.
# Parameters:
#   df = The data frame to be written.
#   fileName = A string representing the name of the file to be written to.
def writeDf(df, fileName):
    outputFile = open(file = fileName, mode = "w")
    outputFile.write(df.to_csv(sep = '\t', header = True))
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
        countFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.counts.txt"
        sectorProbFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.sectorProbs.txt"
        AAProbFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.AAProbs.txt"
        # Last preparations for writing data frames.
        sectorProbDf = countDf.div(countDf.sum(axis = 0), axis = 1) # Probability distribution for each sector
        AAProbDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for each contacting AA
        # Writing frames.
        writeDf(countDf, countFileName)
        writeDf(sectorProbDf, sectorProbFileName)
        writeDf(AAProbDf, AAProbFileName)