import re;

# The area limit for what is considered a contact.
CONTACT_THRESHOLD = 20.0;
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

# Purpose: To fetch info from a given line.
# Parameters:
#   line = The line containing the info.
# Return:
#   A list of parameters fetched from the line.
def fetchLineInfo(line):
    lineList = line.split();  # Splitting by any empty space.
    # Skipping non-informative lines.
    if "SAS" in line or "#" in line or line == "\n" or line == "":
        return False;
    # Collecting important information.
    try:
        # Information about the source particle.
        parName = lineList[1];  # Name of the particle (V0, for example)
        parAA = parName[0];  # One-letter AA code for the particle
        parType = "BB" if parName[-1] == "0" else "SC";  # 0 indicates backbone, 1 indicates sidechain
        parRes = int(lineList[2]);  # Residue number of the particle (Like 1)
        # Information about the contact particle.
        contactName = lineList[5];  # Name of the contact particle (L0, for example)
        contactAA = contactName[0];  # One-letter AA code for the contact particle
        contactType = "BB" if contactName[-1] == "0" else "SC";  # 0 indicates backbone, 1 indicates sidechain
        contactRes = int(lineList[6]);  # Residue of the contact particle (Like 2)
        # General contact information.
        contactArea = float(lineList[8]);  # Area of the contact (>25 indicates a contact)
        sectorInfo = lineList[-1];  # Sector of the contact particle relative to the current particle
        # Calculating the difference between the residue positions.
        resDiff = abs(contactRes - parRes);
    except:
        return False;
    # Skipping lines that are not real contacts or are describing non-real particles.
    if contactArea < CONTACT_THRESHOLD:
        return False;
    if parAA == "X" or contactAA == "X":  # Abnormal amino acids - unknown identity.
        return False;
    # Getting sector information.
    if "=" in sectorInfo:
        sectorNum = int(re.search('(?<=\=).+', lineList[-1]).group(0).rstrip().lstrip());
    else:
        sectorNum = int(sectorInfo.rstrip().lstrip());
    # Skipping lines describing sector 0 contacts.
    if sectorNum == 0:
        return False
    # Returning collected information.
    return (parName, parAA, parType, parRes,
            contactName, contactAA, contactType, contactRes,
            contactArea, sectorNum, resDiff);

# Purpose: To filter a contact set, removing any contacts which
# share a sector with another, keeping the contact with the higher
# area.
# Parameters:
#   contactSet = The set of contacts for a particle.
# Return:
#   filterSet = The filtered set of contacts.
def filterSet(contactSet):
    toRemove = set();
    i = 0;
    while i < len(contactSet):
        contact = contactSet[i];
        otherContacts = [x for x in contactSet if x != contact];
        for otherContact in otherContacts:
            if contact[i_SECTOR_NUM] == otherContact[i_SECTOR_NUM]:
                toPop = contact if contact[i_CONTACT_AREA] < otherContact[i_CONTACT_AREA] else otherContact;
                toRemove.add(toPop);
        i += 1;
    filterSet = [x for x in contactSet if x not in toRemove];
    return filterSet;