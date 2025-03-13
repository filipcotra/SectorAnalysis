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
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"];
# Setting a constant for data frame row names.
ROW_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1", "G0", "G1",
             "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1", "M0", "M1", "N0", "N1",
             "P0", "P1", "Q0", "Q1", "R0", "R1", "S0", "S1", "T0", "T1", "V0", "V1",
             "W0", "W1", "Y0", "Y1"];
# Tracking data frames for each sector. One frame for each sector.
sectorDfs = {1: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             2: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             3: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             4: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             5: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             6: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             7: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             8: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             9: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             10: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             11: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
             12: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES)};

# Purpose: To populate the sector fata frames.
# Parameters:
#   parName = The name of the contacting particle.
#   sectorNum = The sector number of the contact particle.
#   otherContacts = The set of contacts not including the current one.
def populateDict(parName, sectorNum, otherSectors):
    parDf = sectorDfs[sectorNum];
    for otherSector in otherSectors:
        sectorCol = f"S{otherSector[2]}";
        parDf.at[parName, sectorCol] += 1;

# Purpose: To collect stats for occupancy analysis from the line data.
# Parameters:
#   contactSet = The set of contacts for the current particle.
def occupancyAnalysis(contactSet):
    global sectorDfs;
    # Iterating through the contact set.
    for contact in contactSet:
        otherSectors = [x for x in contactSet if x != contact]; # All other contacts
        sourcePar = contact[i_PAR_NAME];
        sourceSector = contact[i_SECTOR_NUM]; # The sector of the current iterable
        sourceDf = sectorDfs[sourceSector];
        for otherSector in otherSectors:
            sectorCol = f"S{otherSector[i_SECTOR_NUM]}";
            sourceDf.at[sourcePar, sectorCol] += 1;

# Purpose: To print the occupancy data.
def printOccupancyOutput():
    global sectorDfs;
    for sector in sectorDfs.keys():
        countDf = sectorDfs[sector];
        countFile = f"sectorOccupancy/S{sector}.occupancy.count.txt";
        probFile = f"sectorOccupancy/S{sector}.occupancy.prob.txt";
        count_outputFile = open(countFile, "w");
        prob_outputFile = open(probFile, "w");
        probDf = countDf.div(countDf.sum(axis = 1), axis = 0); # Along each row
        count_outputFile.write(countDf.to_csv(sep = '\t', header = True));
        prob_outputFile.write(probDf.to_csv(sep = '\t', header = True));
        count_outputFile.close();
        prob_outputFile.close();
