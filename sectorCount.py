import sys
import re
import pandas as pd
from printSectorInfo import printSectorInfo

# Only input should be a file with a list of files.
files = open(sys.argv[1], "r")

# Setting a constant for contact threshold.
CONTACT_THRESHOLD = 25.0
# Setting a constant for data frame column names.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]

# Creating dictionaries to track amino acid objects
bb2bb = {} # Permanent BB2BB
sc2bb = {} # Permanent SC2BB
bb2sc = {} # Permanent BB2SC - inverse of previous
bb_np = {} # Non-permanent BB contacts - with any other particle type
sc_np = {} # Non-permanent SC contacts - with any other particle type

# This function will just populate a given dictionary with a provided value.
# Parameters:
#   parKey = A string representing the key for the current particle.
#   contactKey = A string representing the key for the contacting particle.
#   sectoryKey = A string representing the sector key (S1, S2, etc.).
#   dict = the dictionary to be filled.
def populateDict(parKey, contactKey, sectorKey, dict):
    if parKey in dict.keys():
        # If the contacted particle key has been encountered before, add to it.
        if contactKey in dict[parKey].index:
            dict[parKey].at[contactKey, sectorKey] += 1
        else:
            dict[parKey].loc[contactKey] = 0 # Adding a new row where everything is equal to 0
            dict[parKey].at[contactKey, sectorKey] += 1
    # If the parKey has not been encountered, everything has to be initialized.
    else:
        dict[parKey] = pd.DataFrame(columns = COL_NAMES)
        dict[parKey].loc[contactKey] = 0 # Adding a new row where everything is equal to 0
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
            parType = "bb" if parName[-1] == "0" else "sc" # 0 indicates backbone, 1 indicates sidechain
            resNum = int(lineList[2]) # Residue number of the particle (Like 1)
            # Information about the contact particle.
            contactName = lineList[5] # Name of the contact particle (L0, for example)
            contactAA = contactName[0] # One-letter AA code for the contact particle
            contactType = "bb" if contactName[-1] == "0" else "sc" # 0 indicates backbone, 1 indicates sidechain
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
        parKey = f"{parAA}.{parType}"
        contactKey = f"{contactAA}.{contactType}"
        sectorKey = f"S{sectorNum}"
        # Populating contacts
        if parType == "bb":
            if (resDiff == 1 and contactType == "bb"): # Indicating a BB2BB contact
                populateDict(parKey, contactKey, sectorKey, bb2bb)
            elif (resDiff == 0 and contactType == "sc"): # Indicating a BB2SC contact
                populateDict(parKey, contactKey, sectorKey, bb2sc)
            else: # Indicates a non-permanent BB contact
                populateDict(parKey, contactKey, sectorKey, bb_np)
        else: # Indicating an SC particle
            if (resDiff == 0 and contactType == "bb"): # Indicating a SC2BB contact
                populateDict(parKey, contactKey, sectorKey, sc2bb)
            else: # Indicates a non-permanent SC contact
                populateDict(parKey, contactKey, sectorKey, sc_np)
    currFile.close()
# Writing results to files.
printSectorInfo(bb2bb, "bb2bb")
printSectorInfo(bb2sc, "bb2sc")
printSectorInfo(sc2bb, "sc2bb")
printSectorInfo(bb_np, "bb_np")
printSectorInfo(sc_np, "sc_np")