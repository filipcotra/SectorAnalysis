import pandas as pd;

# Setting a constant for data frame column names.
COL_NAMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"];
# Tracking the dataframes made.
contactDict = {};

# Purpose: Populate a given dictionary with a provided value.
# Parameters:
#   parName = A string representing the name of the current particle.
#   contactNum = The number of contacts the current particle has.
def populateDict(parName, contactNum):
    global contactDict;
    sectorCol = f"{contactNum}";
    # If the amino-acid has not been encountered before
    if parName not in contactDict.keys():
        contactDict[parName] = pd.DataFrame(columns = COL_NAMES);
        contactDict[parName].loc["Count"] = 0;
    # Update the appropriate data frame, iterating to indicate the number
    # of contacts.
    contactDict[parName].at["Count", sectorCol] += 1;

# Tracking the current particle and the number of contacts it has.
currPar = None;
currCount = 0;
# Purpose: To track appropriate capacity stats in data frames.
# Parameters:
#   parName = The name of the given particle (XN, where X is the AA and N is 0 or 1).
#   parRes = The residue number of the given particle.
#   contactName = The name of the contact particle (XN, where X is the AA and N is 0 or 1).
#   contactAA = The one-letter AA code for the contact particle.
#   sectorNum = The sector of the current contact.
def capacityAnalysis(parName, parRes):
    global currPar, currCount;
    parKey = (parName, parRes);
    if currPar is None:
        currPar = parKey;
    # If the new particle is not the same as the current one, that means we've moved on.
    # Should store current particle information and reassign the new particle.
    if parKey != currPar:
        currName = currPar[0];
        populateDict(currName, currCount);
        currCount = 0;
        currPar = parKey;
    currCount += 1;

# Purpose: To print the output.
def printParticleOutput():
    global contactDict;
    for resName in contactDict.keys():
        count_fileName = f"particleCapacity/{resName}.contactCapacity.count.txt"
        prob_fileName = f"particleCapacity/{resName}.contactCapacity.prob.txt"
        countDf = contactDict[resName]
        probDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for each contacting AA
        probDf = probDf.rename(index={"Count": "Prob"})
        # Writing to files.
        countFile = open(count_fileName, "w")
        probFile = open(prob_fileName, "w")
        countFile.write(countDf.to_csv(sep='\t', header=True))
        probFile.write(probDf.to_csv(sep='\t', header=True))


