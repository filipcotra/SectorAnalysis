import sys
import re
import pandas as pd

# Setting a constant for contact threshold.
CONTACT_THRESHOLD = 25.0
# Only input should be a file with a list of files.
files = open(sys.argv[1], "r")
# Setting a constant for data frame column names.
COL_NAMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
# Tracking the dataframes made.
contactDict = {}

# This function will just populate a given dictionary with a provided value.
# Parameters:
#   parName = A string representing the name of the current particle.
#   contactNum = A string representing the number of contacts the current particle has.
def populateDict(parName, contactNum):
    contactKey = f"{contactNum}"
    # If the amino-acid has been encountered before
    if parName in contactDict.keys():
        contactDict[parName].at["Count", contactKey] += 1
    # If the parName has not been encountered, everything has to be initialized.
    else:
        contactDict[parName] = pd.DataFrame(columns = COL_NAMES)
        contactDict[parName].loc["Count"] = 0
        contactDict[parName].at["Count", contactKey] += 1

# Tracking current particle and contact count.
currPar = None;
currCount = 0;
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
            if currPar is None:
                currPar = parName
            parAA = parName[0] # One-letter AA code for the particle
            # Information about the contact particle.
            contactName = lineList[5] # Name of the contact particle (L0, for example)
            contactAA = contactName[0] # One-letter AA code for the contact particle
            # General contact information.
            contactArea = float(lineList[8]) # Area of the contact (>25 indicates a contact)
            sectorInfo = lineList[-1] # Sector of the contact particle relative to the current particle
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
        # If the new particle is not the same as the current one, that means we've moved on.
        # Should store current particle information and reassign the new particle.
        if parName != currPar:
            populateDict(currPar, currCount)
            currCount = 0
            currPar = parName
        currCount += 1
    currFile.close()

# Printing output.
for resName in contactDict.keys():
    count_fileName = f"{resName}.contactCapacity.count.txt"
    prob_fileName = f"{resName}.contactCapacity.prob.txt"
    countDf = contactDict[resName]
    probDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for each contacting AA
    probDf = probDf.rename(index={"Count": "Prob"})
    # Writing to files.
    countFile = open(count_fileName, "w")
    probFile = open(prob_fileName, "w")
    countFile.write(countDf.to_csv(sep='\t', header=True))
    probFile.write(probDf.to_csv(sep='\t', header=True))


