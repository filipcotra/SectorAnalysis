import pandas as pd;

# Setting a constant for data frame column names.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]
ROW_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1", "G0", "G1",
             "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1", "M0", "M1", "N0", "N1",
             "P0", "P1", "Q0", "Q1", "R0", "R1", "S0", "S1", "T0", "T1", "V0", "V1",
             "W0", "W1", "Y0", "Y1"];
# Creating dictionaries to track amino acid objects
bb2bb = {}; # Permanent BB2BB
sc2bb = {}; # Permanent SC2BB
bb2sc = {}; # Permanent BB2SC - inverse of previous
bb_np = {}; # Non-permanent BB contacts - with any other particle type
sc_np = {}; # Non-permanent SC contacts - with any other particle type

# Purpose: Simple function to write to file.
# Parameters:
#   df = The data frame to be written.
#   fileName = A string representing the name of the file to be written to.
def writeDf(df, fileName):
    outputFile = open(file = fileName, mode = "w");
    outputFile.write(df.to_csv(sep = '\t', header = True));
    outputFile.close();

# Purpose: Populate a given dictionary with a provided value.
# Parameters:
#   parKey = A string representing the key for the current particle.
#   contactKey = A string representing the key for the contacting particle.
#   sectoryKey = A string representing the sector key (S1, S2, etc.).
#   dict = the dictionary to be filled.
def populateDict(parKey, contactKey, sectorKey, dict):
    if parKey in dict.keys():
        dict[parKey].at[contactKey, sectorKey] += 1
    # If the parKey has not been encountered, everything has to be initialized.
    else:
        dict[parKey] = pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES)
        dict[parKey].at[contactKey, sectorKey] += 1

# Purpose: To collect data relevant to sector counts.
# Parameters:
#   parName = The name of the given particle.
#   parType = The type of the given particle.
#   contactName = The name of the contact particle.
#   contactType = The type of the contact particle.
#   resDiff = The difference between the residue numbers.
#   sectorNum = The sector that the contact is within.
def countAnalysis(parName, parType, contactName, contactType, resDiff, sectorNum):
    global bb2bb, bb2sc, sc2bb, bb_np, sc_np;
    # This will be used as the key for the dictionaries
    sectorKey = f"S{sectorNum}"
    # Populating contacts
    if parType == "BB":
        if (resDiff == 1 and contactType == "BB"): # Indicating a BB2BB contact
            populateDict(parName, contactName, sectorKey, bb2bb)
        elif (resDiff == 0 and contactType == "SC"): # Indicating a BB2SC contact
            populateDict(parName, contactName, sectorKey, bb2sc)
        else: # Indicates a non-permanent BB contact
            populateDict(parName, contactName, sectorKey, bb_np)
    else: # Indicating an SC particle
        if (resDiff == 0 and contactType == "BB"): # Indicating a SC2BB contact
            populateDict(parName, contactName, sectorKey, sc2bb)
        else: # Indicates a non-permanent SC contact
            populateDict(parName, contactName, sectorKey, sc_np)

# Purpose: Take a dictionary and print out all of its information
# into separate files.
# Parameters:
#   dict = The dictionary containing dataframe objects.
#   dictType = The name of the dictionary (bb2bb, bb2sc, etc.).
def printSectorInfo(dict, dictName):
    for parKey in dict.keys():
        countDf = dict[parKey]; # Dataframes stored in the dictionary have count values
        # Making a file names.
        countFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.counts.txt";
        sectorProbFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.sectorProbs.txt";
        AAProbFileName = f"sectorCount/{dictName}/{parKey}.{dictName}.AAProbs.txt";
        # Last preparations for writing data frames.
        sectorProbDf = countDf.div(countDf.sum(axis = 0), axis = 1); # Probability distribution for each sector
        AAProbDf = countDf.div(countDf.sum(axis = 1), axis = 0); # Probability distribution for each contacting AA
        # Replacing any NaNs with 0s.
        sectorProbDf = sectorProbDf.fillna(0);
        AAProbDf = AAProbDf.fillna(0);
        # Writing frames.
        writeDf(countDf, countFileName);
        writeDf(sectorProbDf, sectorProbFileName);
        writeDf(AAProbDf, AAProbFileName);

# Purpose: Call the printSectorInfo function with the relevant
# dictionaries.
def printCountOutput():
    global bb2bb, bb2sc, sc2bb, bb_np, sc_np;
    # Writing results to files.
    printSectorInfo(bb2bb, "bb2bb");
    printSectorInfo(bb2sc, "bb2sc");
    printSectorInfo(sc2bb, "sc2bb");
    printSectorInfo(bb_np, "bb_np");
    printSectorInfo(sc_np, "sc_np");