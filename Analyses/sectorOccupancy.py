import pandas as pd;

# Setting a constant for data frame column names.
COL_NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"];
# Setting a constant for data frame row names.
ROW_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1", "G0", "G1",
             "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1", "M0", "M1", "N0", "N1",
             "P0", "P1", "Q0", "Q1", "R0", "R1", "S0", "S1", "T0", "T1", "V0", "V1",
             "W0", "W1", "Y0", "Y1"];
# Tracking data frames for each sector. One frame for each sector.
sectorDFs = {1: pd.DataFrame(0, columns = COL_NAMES, index = ROW_NAMES),
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

# Tracking current particle and contact.
currPar = None;
currContacts = [];
# Purpose: To collect stats for occupancy analysis from the line data.
# Parameters:
#   parName = The name of the given particle.
#   parRes = The residue number for the given particle.
#   contactName = The name of the contact particle.
#   contactRes = The residue number for the contact particle.
#   sectorNum = The sector number for the given contact.
def occupancyAnalysis(parName, parRes, contactName, contactRes, sectorNum):
    global sectorDFs, currPar, currContacts;
    parKey = (parName, parRes);
    # If the new particle is not the same as the current one, that means we've moved on.
    # Should store current particle information and reassign the new particle.
    if parKey != currPar:
        # Populating.
        for contact in currContacts:
            otherSectors = [x for x in currContacts if x != contact]; # All other contacts
            sourcePar = contact[0];
            sourceSector = contact[2]; # The sector of the current iterable
            sourceDf = sectorDFs[sourceSector];
            for otherSector in otherSectors:
                sectorCol = f"S{otherSector[2]}";
                sourceDf.at[sourcePar, sectorCol] += 1;
        currContacts = [];
        currPar = parKey;
    # Regardless, we should populate the current contacts.
    currContacts.append((contactName, contactRes, sectorNum)); # Ex: (V0, 1, 3).

# Purpose: To print the occupancy data.
def printOccupancyOutput():
    global sectorDFs;
    for sector in sectorDFs.keys():
        countDf = sectorDFs[sector];
        countFile = f"sectorOccupancy/S{sector}.occupancy.count.txt";
        probFile = f"sectorOccupancy/S{sector}.occupancy.prob.txt";
        count_outputFile = open(countFile, "w");
        prob_outputFile = open(probFile, "w");
        probDf = countDf.div(countDf.sum(axis = 1), axis = 0); # Along each row
        count_outputFile.write(countDf.to_csv(sep = '\t', header = True));
        prob_outputFile.write(probDf.to_csv(sep = '\t', header = True));
        count_outputFile.close();
        prob_outputFile.close();
