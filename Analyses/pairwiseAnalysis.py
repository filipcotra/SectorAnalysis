import pandas as pd;
import numpy as np;

# Defining constants for navigating contact sets.
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
i_CONTACT_DISTANCE = 11;
# Defining constants for navigating edges.
i_EDGE_NAME_A = 0;
i_EDGE_RES_A = 1;
i_EDGE_SECTOR_A = 2;
i_EDGE_NAME_B = 3;
i_EDGE_RES_B = 4;
i_EDGE_SECTOR_B = 5;
# Defining particle names.
PAR_NAMES = ["A0", "A1", "C0", "C1", "D0", "D1", "E0", "E1", "F0", "F1",
             "G0", "G1", "H0", "H1", "I0", "I1", "K0", "K1", "L0", "L1",
             "M0", "M1", "N0", "N1", "P0", "P1", "Q0", "Q1", "R0", "R1",
             "S0", "S1", "T0", "T1", "V0", "V1", "W0", "W1", "Y0", "Y1"];
SECTORS = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"];

# Making data frames to hold the individual stats for non-permanent contacts.
# None of the data frames will have indexes - they must be added as they are
# identified. No redundancies will be allowed.
# ------------------- One stat -------------------
df_parSec_NP = pd.DataFrame(columns = ["Count"]); # Particle/Sector to Particle/Sector pair incidence
df_parCap_NP = pd.DataFrame(data = np.zeros((len(PAR_NAMES), 13)),
                            index = PAR_NAMES,
                            columns = [f"{x}" for x in range(13)]); # The number of contacts a particle has
df_sharedCon_NP = pd.DataFrame(columns = ["Count"]); # The number of shared contacts N between the two particles of a contact
df_seqSep_NP = pd.DataFrame(columns = ["Count"]); # The sequence separation between two contacting particles
df_conSep_NP = pd.DataFrame(columns = ["Count"]); # The number of contacts separating two particles, ignoring the current edge.
# ------------------- Five stats -------------------
# Particle/sector pairwise x shared contacts x sequence separation x particle capacity x contact separation.
df_all_NP = pd.DataFrame(columns = ["Particle.Sector.Pair", "Shared.Contacts", "Sequence.Separation",
                                    "Particle.Capacity", "Contact.Separation", "Count"]);

# Making data frames to hold the individual stats for permanent contacts.
# These will not be separated for permanent contact type, as it should
# be fairly obvious from the resulting data (which the model should account for).
df_P = pd.DataFrame(columns = ["Particle.Sector.Pair", "Sequence.Difference", "Count"]); # Particle pairs may influence sector pairs.

# Purpose: To analyze stats pertaining to contact pairs and
# particle properties, particle types, sector interactions,
# sequence separation, and number of shared contacts.
# Parameters:
#   edgeSet = The set of edges describing the structure.
#   contactSet = A dictionary containing the contact set for
#                each particle in the structure.
def analyze(edgeSet, contactSet):
    for edge in edgeSet:
        resDiff = abs(edge[i_EDGE_RES_A] - edge[i_EDGE_RES_B]);
        parA_isBB = edge[i_EDGE_NAME_A][-1] == "0";
        parB_isBB = edge[i_EDGE_NAME_B][-1] == "0";
        # Only happens with R group bonds or peptide bonds.
        if resDiff == 0 or (resDiff == 1 and parA_isBB and parB_isBB):
            analyzeP(edge);
        else:
            analyzeNP(edge, contactSet);

# Purpose: To populate data frames with non-permanent contact
# information.
# Parameters:
#   edge = The current edge being analyzed.
#   edgeSet = The set of edges describing the structure.
#   contactSet = The set of contact sets in the structure.
def analyzeNP(edge, contactSet):
    # Particle A information.
    A_parName = edge[i_EDGE_NAME_A];
    A_parRes = edge[i_EDGE_RES_A];
    A_sector = edge[i_EDGE_SECTOR_A];
    A_par = (A_parName, A_parRes);
    A_contactSet = contactSet[A_par];
    A_numContacts = len(A_contactSet);
    # Particle B information.
    B_parName = edge[i_EDGE_NAME_B];
    B_parRes = edge[i_EDGE_RES_B];
    B_sector = edge[i_EDGE_SECTOR_B];
    B_par = (B_parName, B_parRes);
    B_contactSet = contactSet[B_par];
    B_numContacts = len(B_contactSet);
    # Information pertaining to both particles.
    seqDiff = A_parRes - B_parRes;
    seqSep = abs(seqDiff);
    shared = getNumShared(A_contactSet, B_contactSet);
    numShared = len(shared);
    contactSep = getContactSep(A_par, B_par, seqSep, A_sector, contactSet);
    # Setting some variables to define key orders. If done
    # properly, this should remove the need to check for
    # reverse keys.
    if A_parName[-1] == B_parName[-1]: # Both particles are of the same type.
        parName_1 = A_parName if A_parName[0] <= B_parName[0] else B_parName;
    else: # Should always be 0 before 1.
        parName_1 = A_parName if A_parName[-1] <= B_parName[-1] else B_parName;
    parName_2 = B_parName if parName_1 is A_parName else A_parName;
    sector_1 = A_sector if parName_1 is A_parName else B_sector;
    sector_2 = B_sector if sector_1 is A_sector else A_sector;
    parCap_1 = A_numContacts if parName_1 is A_parName else B_numContacts;
    parCap_2 = B_numContacts if parCap_1 is A_numContacts else A_numContacts;
    # Collecting stats.
    # ----- One stat -----
    # Building keys.
    key_parSec = f"{parName_1}/S{sector_1};{parName_2}/S{sector_2}";
    key_parCap = f"{parCap_1};{parCap_2}";
    key_sharedCon = f"{numShared}";
    key_seqSep = f"{seqSep}";
    key_conSep = f"{contactSep}";
    # Counting.
    populateDf(rowKey = key_parSec, df = df_parSec_NP);
    populateDf(rowKey = parName_1, df = df_parCap_NP, colKey = f"{parCap_1}"); # Particle A capacity.
    populateDf(rowKey = parName_2, df = df_parCap_NP, colKey = f"{parCap_2}"); # Particle B capacity.
    populateDf(rowKey = key_sharedCon, df = df_sharedCon_NP);
    populateDf(rowKey = key_seqSep, df = df_seqSep_NP);
    populateDf(rowKey = key_conSep, df = df_conSep_NP);
    # ----- Five stat keys -----
    # Building key.
    key_all = "|".join([key_parSec, key_parCap, key_sharedCon, key_seqSep, key_conSep]);
    # Counting.
    if key_all in df_all_NP.index:
        df_all_NP.loc[key_all, "Count"] += 1;
    else:
        df_all_NP.loc[key_all] = [key_parSec,  key_parCap, key_sharedCon, key_seqSep, key_conSep, 1];

# Purpose: To populate data frames with permanent contact
# information.
# Parameters:
#   edge = The current edge being analyzed.
def analyzeP(edge):
    # Particle A information.
    A_parName = edge[i_EDGE_NAME_A];
    A_parRes = edge[i_EDGE_RES_A];
    A_sector = edge[i_EDGE_SECTOR_A];
    # Particle B information.
    B_parName = edge[i_EDGE_NAME_B];
    B_parRes = edge[i_EDGE_RES_B];
    B_sector = edge[i_EDGE_SECTOR_B];
    # Information pertaining to both particles.
    seqDiff = A_parRes - B_parRes;
    # Setting some variables to define key orders. If done
    # properly, this should remove the need to check for
    # reverse keys.
    if A_parName[-1] == B_parName[-1]:  # Both particles are of the same type.
        parName_1 = A_parName if A_parName[0] <= B_parName[0] else B_parName;
    else:  # Should always be 0 before 1.
        parName_1 = A_parName if A_parName[-1] <= B_parName[-1] else B_parName;
    parName_2 = B_parName if parName_1 is A_parName else A_parName;
    sector_1 = A_sector if parName_1 is A_parName else B_sector;
    sector_2 = B_sector if sector_1 is A_sector else A_sector;
    # Collecting stats.
    # Building keys.
    key_parSec = f"{parName_1}/S{sector_1};{parName_2}/{sector_2}";
    key_seqDiff = f"{seqDiff}";
    # Building keys.
    key_all = "|".join([key_parSec, key_seqDiff]);
    # Counting.
    if key_all in df_P.index:
        df_P.loc[key_all, "Count"] += 1;
    else:
        df_P.loc[key_all] = [key_parSec, key_seqDiff, 1];

# Purpose: To find the contact separation between two particles.
# Parameters:
#   A_par = Particle A.
#   B_par = Particle B.
#   ignoreSector_A = The sector of particle A known to contain particle B,
#                    which should thereby be ignored.
#   seqSep = The sequence separation between the two particles.
#   contactSet = The set of contact sets in the structure.
def getContactSep(A_par, B_par, ignoreSector_A, seqSep, contactSet):
    sep = 0;
    maxSep = seqSep - 1;
    checkedPars = set();
    parsToCheck = {A_par};
    # Every particle is connected to every other particle, so
    # eventually a number will be reported. The minimum number
    # can actually be found by using the residue difference,
    # which is equal to 1 greater than the number of separating
    # edges given a linear chain.
    while sep < maxSep: # If sep ever equals seqSep, return sep.
        contactPars = set();
        for par in parsToCheck: # Iterating through all the pars that need to be checked - contacts of the previous par.
            if par is A_par:
                parContacts = [(x[i_CONTACT_NAME], x[i_CONTACT_RES]) for x in contactSet[par] if x[i_SECTOR_NUM] != ignoreSector_A];
            else:
                parContacts = [(x[i_CONTACT_NAME], x[i_CONTACT_RES]) for x in contactSet[par]];
            filteredContacts = [x for x in parContacts if x not in checkedPars]; # All contacting pars that have not been checked already.
            if B_par in filteredContacts: # If par B is in the contact set, return the current sep value.
                return sep;
            contactPars.update(filteredContacts);
        # All of the particles in parsToCheck have been checked.
        checkedPars.update(parsToCheck);
        # Updating parsToCheck with the new list of contacting particles.
        parsToCheck = contactPars;
        # Indicating that one edge now separates particles A and B.
        sep += 1;
        del contactPars;
    return sep;

# Purpose: To find the set of shared contacts and return their
# amount.
# Parameters:
#   A_contactSet = Contact set A.
#   B_contactSet = Contact set B.
# Return:
#   shared = A list of particles in contact with both A and B.
#   total = A list of particles in contact with A or B.
def getNumShared(A_contactSet, B_contactSet):
    A_contacts = set([(x[i_CONTACT_NAME], x[i_CONTACT_RES]) for x in A_contactSet]);
    B_contacts = set([(x[i_CONTACT_NAME], x[i_CONTACT_RES]) for x in B_contactSet]);
    shared = list(A_contacts & B_contacts);
    return shared;

# Purpose: To populate a given data frame.
# Parameters:
#   rowKey = The key for the row to populate.
#   df = The data frame to populate.
#   colKey = The key for the column to populate. Default is "Count".
def populateDf(rowKey, df, colKey = "Count"):
    if rowKey in df.index:
        df.loc[rowKey, colKey] += 1;
    else:
        df.loc[rowKey, colKey] = 1;

# Purpose: To write stat dataframes to files.
# Parameters:
#   outputPath = The path to the output file.
def export(outputPath):
    # Non-permanent stats.
    # ------------------- One stat -------------------
    df_parSec_NP.to_csv(f"{outputPath}/NP_stats_particleSectorPairs.tsv", sep = "\t", index = True);
    df_parCap_NP.to_csv(f"{outputPath}/NP_stats_particleCapacity.tsv", sep = "\t", index = True);
    df_sharedCon_NP.to_csv(f"{outputPath}/NP_stats_sharedContacts.tsv", sep = "\t", index = True);
    df_seqSep_NP.to_csv(f"{outputPath}/NP_stats_sequenceSeparation.tsv", sep = "\t", index = True);
    df_conSep_NP.to_csv(f"{outputPath}/NP_stats_contactSeparation.tsv", sep = "\t", index = True);
    # ------------------- Five stats -------------------
    df_all_NP.to_csv(f"{outputPath}/NP_stats_all.tsv", sep = "\t", index = True);
    # Permanent stats
    df_P.to_csv(f"{outputPath}/P_stats.tsv", sep = "\t", index = True);