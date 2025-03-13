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
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"]
ROW_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1", "G0", "G1",
             "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1", "M0", "M1", "N0", "N1",
             "P0", "P1", "Q0", "Q1", "R0", "R1", "S0", "S1", "T0", "T1", "V0", "V1",
             "W0", "W1", "Y0", "Y1"];
# Creating dictionaries to track amino acid objects
sectorCounts_NP = {};
sectorCounts_peptideBonds_minus1 = {};
sectorCounts_peptideBonds_plus1 = {};
sectorCounts_rGroupBonds = {};


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
#   sectorKey = A string representing the sector key (S1, S2, etc.).
#   dict = the dictionary to be filled.
def populateDict(parKey, contactKey, sectorKey, dict):
    if parKey not in dict.keys():
        dict[parKey] = pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES);
    dict[parKey].at[contactKey, sectorKey] += 1

# Purpose: To collect data relevant to sector counts for non-permanent sector counts.
# Parameters:
#   contactSet = The set of contacts for the current particle.
def countAnalysis(contactSet):
    global sectorCounts_NP;
    # Iterate through the contact set.
    for contact in contactSet:
        sectorKey = f"S{contact[i_SECTOR_NUM]}";
        parName = contact[i_PAR_NAME];
        parType = contact[i_PAR_TYPE];
        parRes = contact[i_PAR_RES];
        contactName = contact[i_CONTACT_NAME];
        contactType = contact[i_CONTACT_TYPE];
        contactRes = contact[i_CONTACT_RES];
        resDiff = contact[i_RES_DIFF];
        # Populating contacts
        if parType == "BB":
            if (resDiff == 1 and contactType == "BB"): # Indicating a permanent BB2BB contact - skip.
                if parRes > contactRes:
                    populateDict(parName, contactName, sectorKey, sectorCounts_peptideBonds_minus1);
                else:
                    populateDict(parName, contactName, sectorKey, sectorCounts_peptideBonds_plus1);
            elif (resDiff == 0 and contactType == "SC"): # Indicating a permanent BB2SC contact - skip.
                populateDict(parName, contactName, sectorKey, sectorCounts_rGroupBonds);
            else: # Indicates a non-permanent BB contact
                populateDict(parName, contactName, sectorKey, sectorCounts_NP);
        else: # Indicating an SC particle
            if (resDiff == 0 and contactType == "BB"): # Indicating a permanent SC2BB contact - skip.
                populateDict(parName, contactName, sectorKey, sectorCounts_rGroupBonds);
            else: # Indicates a non-permanent SC contact
                populateDict(parName, contactName, sectorKey, sectorCounts_NP);

# Purpose: Print dictionary information to files.
# Parameters:
#   dict = The dictionary to print out.
def printSectorInfo(dict):
    global sectorCounts_NP, sectorCounts_peptideBonds_minus1, sectorCounts_peptideBonds_plus1, sectorCounts_rGroupBonds;
    fileSig = None;
    if dict is sectorCounts_NP:
        fileSig = "NP";
    elif dict is sectorCounts_peptideBonds_minus1:
        fileSig = "BB2BB_minus1";
    elif dict is sectorCounts_peptideBonds_plus1:
        fileSig = "BB2BB_plus1";
    elif dict is sectorCounts_rGroupBonds:
        fileSig = "rGroup"
    for parKey in dict.keys():
        countDf = dict[parKey]; # Dataframes stored in the dictionary have count values
        sectorProbDf = countDf.div(countDf.sum(axis = 0), axis = 1); # Probability distribution for each sector
        AAProbDf = countDf.div(countDf.sum(axis = 1), axis = 0); # Probability distribution for each AA
        # Replacing any NaNs with 0s.
        sectorProbDf = sectorProbDf.fillna(0);
        AAProbDf = AAProbDf.fillna(0);
        # Making a file names.
        countFile = f"sectorCount/{parKey}.{fileSig}.sectorCounts.txt";
        sectorProbFile = f"sectorCount/{parKey}.{fileSig}.sectorProbs.txt";
        AAProbFile = f"sectorCount/{parKey}.{fileSig}.AAProbs.txt";
        # Opening files.
        count_outputFile = open(countFile, "w");
        sectorProb_outputFile = open(sectorProbFile, "w");
        AAProb_outputFile = open(AAProbFile, "w");
        # Writing to files.
        count_outputFile.write(countDf.to_csv(sep = '\t', header = True));
        sectorProb_outputFile.write(sectorProbDf.to_csv(sep = '\t', header = True));
        AAProb_outputFile.write(AAProbDf.to_csv(sep = '\t', header = True));

# Purpose: Call the printSectorInfo function with the relevant
# dictionaries.
def printCountOutput():
    printSectorInfo(sectorCounts_NP);
    printSectorInfo(sectorCounts_peptideBonds_minus1);
    printSectorInfo(sectorCounts_peptideBonds_plus1);
    printSectorInfo(sectorCounts_rGroupBonds);