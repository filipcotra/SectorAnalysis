import sys
import re
import pandas as pd

# Only input should be a file with a list of files.
files = open(sys.argv[1], "r")

# Setting a constant for contact threshold.
CONTACT_THRESHOLD = 25.0
# Setting a constant for data frame column names.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]
ROW_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1", "G0", "G1",
             "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1", "M0", "M1", "N0", "N1",
             "P0", "P1", "Q0", "Q1", "R0", "R1", "S0", "S1", "T0", "T1", "V0", "V1",
             "W0", "W1", "Y0", "Y1"]
# Creating dictionaries to track amino acid objects
bb2bb = {} # Permanent BB2BB
sc2bb = {} # Permanent SC2BB
bb2sc = {} # Permanent BB2SC - inverse of previous
bb_np = {} # Non-permanent BB contacts - with any other particle type
sc_np = {} # Non-permanent SC contacts - with any other particle type

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
        # Replacing any NaNs with 0s.
        sectorProbDf = sectorProbDf.fillna(0);
        AAProbDf = AAProbDf.fillna(0);
        # Writing frames.
        writeDf(countDf, countFileName)
        writeDf(sectorProbDf, sectorProbFileName)
        writeDf(AAProbDf, AAProbFileName)

# This function will just populate a given dictionary with a provided value.
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

# Looping through the files, opening them.
for fileName in files:
    currFile = open(fileName.rstrip(), "r")
    # In this file, identifying the sectors for non-backbone contacts and
    # adding the respective amino acids to the sector values.
    for line in currFile:
        lineList = line.split() # Splitting by any empty space.
        # Skipping non-informative lines.
        if "SAS" in line or "#" in line or line == "\n" or line == "":
            continue
        # Collecting important information.
        try:
            # Information about the source particle.
            parName = lineList[1] # Name of the particle (V0, for example)
            parAA = parName[0] # One-letter AA code for the particle
            parType = "BB" if parName[-1] == "0" else "SC" # 0 indicates backbone, 1 indicates sidechain
            resNum = int(lineList[2]) # Residue number of the particle (Like 1)
            # Information about the contact particle.
            contactName = lineList[5] # Name of the contact particle (L0, for example)
            contactAA = contactName[0] # One-letter AA code for the contact particle
            contactType = "BB" if contactName[-1] == "0" else "SC" # 0 indicates backbone, 1 indicates sidechain
            contactResNum = int(lineList[6]) # Residue of the contact particle (Like 2)
            # General contact information.
            contactArea = float(lineList[8]) # Area of the contact (>25 indicates a contact)
            sectorInfo = lineList[-1] # Sector of the contact particle relative to the current particle
            # Calculating the difference between the residue positions.
            resDiff = abs(contactResNum - resNum)
        except:
            continue
        if parAA == "X" or contactAA == "X": # Abnormal amino acids - unknown identity.
            continue
        # Getting the sector information.
        if "=" in sectorInfo:
            sectorNum = int(re.search('(?<=\=).+', lineList[-1]).group(0).rstrip().lstrip())
        else:
            sectorNum = int(sectorInfo.rstrip().lstrip())
        # Tracking any s = 0 values. If found, just move on.
        if sectorNum == 0:
            continue
        # Tracking contact information.
        if contactArea < CONTACT_THRESHOLD:
            continue
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
    currFile.close()
# Writing results to files.
printSectorInfo(bb2bb, "bb2bb")
printSectorInfo(bb2sc, "bb2sc")
printSectorInfo(sc2bb, "sc2bb")
printSectorInfo(bb_np, "bb_np")
printSectorInfo(sc_np, "sc_np")