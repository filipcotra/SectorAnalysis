import pandas as pd;

# Defining constants for fetchInfo returns.
i_PAR_NAME = 0;
i_PAR_AA = 1;
i_PAR_TYPE = 2;
i_PAR_RES = 3;
i_CONTACT_NAME = 4;
i_CONTACT_AA = 5;
i_CONTACT_TYPE = 6;
i_CONTACT_RES = 7;
i_CONTACT_AREA = 8;
i_SECTOR_NUM = 9;
i_RES_DIFF = 10;
# Setting a constant for data frame column names.
COL_NAMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"];
# Tracking the dataframes made.
contactDict = {};

# Purpose: Populate a given dictionary with a provided value.
# Parameters:
#   parName = A string representing the name of the current particle.
#   contactNum = The number of contacts that the particle has.
def populateDict(parName, contactNum):
    global contactDict;
    contactCol = f"{contactNum}";
    # If the amino-acid has not been encountered before
    if parName not in contactDict.keys():
        contactDict[parName] = pd.DataFrame(columns = COL_NAMES);
        contactDict[parName].loc["Count"] = 0;
    # Update the appropriate data frame, iterating to indicate the number
    # of contacts.
    contactDict[parName].at["Count", contactCol] += 1;

# Purpose: To track appropriate capacity stats in data frames.
# Parameters:
#   currPar = The current particle.
#   contactSet = The set of contacts for the current particle.
def capacityAnalysis(currPar, contactSet):
    parName = currPar[i_PAR_NAME];
    # Iterating through the contactSet to populate the data frames.
    contactNum = len(contactSet);
    populateDict(parName, contactNum);

# Purpose: To print the output.
def printParticleOutput():
    global contactDict;
    for resName in contactDict.keys():
        count_fileName = f"particleCapacity/{resName}.contactCapacity.count.txt"
        prob_fileName = f"particleCapacity/{resName}.contactCapacity.prob.txt"
        countDf = contactDict[resName]
        probDf = countDf.div(countDf.sum(axis = 1), axis = 0) # Probability distribution for number of contacts.
        probDf = probDf.rename(index = {"Count": "Prob"})
        # Writing to files.
        countFile = open(count_fileName, "w")
        probFile = open(prob_fileName, "w")
        countFile.write(countDf.to_csv(sep = '\t', header=True))
        probFile.write(probDf.to_csv(sep = '\t', header=True))


