import numpy as np;

# Index constants for keys.
i_KEY_NAME = 0;
i_KEY_RES = 1;
# Defining constants for contactSet constants.
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

# Purpose: Populating the contents of a given file's
# contact sets into an edge set and returning it.
# Parameters:
#   contactSets = A dictionary with keys representing
#                 particles and values representing
#                 their contact sets as arrays.
# Return:
#   edgeSet = The set of edges in the given contact set.
def getEdgeSet(contactSets):
    edgeSet = set();
    # Iterating through the particles.
    for parKey in contactSets.keys():
        parContacts = contactSets[parKey];
        parName = parKey[i_KEY_NAME];
        parRes = parKey[i_KEY_RES];
        # Iterating through the contacts, making the edges.
        for contact in parContacts:
            parSector = contact[i_SECTOR_NUM];
            # Collecting information about the contact.
            contactName = contact[i_CONTACT_NAME];
            contactRes = contact[i_CONTACT_RES];
            contactKey = (contactName, contactRes);
            # Now that we have the contact key, we can
            # find the sector corresponding to this
            # contact.
            corrSet = contactSets[contactKey];
            corrContact = [x for x in corrSet if x[i_CONTACT_NAME] == parName and x[i_CONTACT_RES] == parRes];
            try:
                contactSector = corrContact[0][i_SECTOR_NUM];
                # Now that we have the corresponding sector,
                # we can make the edge.
                edge = (parName, parRes, parSector, contactName, contactRes, contactSector);
                reverseEdge = (contactName, contactRes, contactSector, parName, parRes, parSector);
                # Making sure not to add redundant information.
                if reverseEdge not in edgeSet:
                    edgeSet.add(edge);
            # No corresponding contact exists, so this edge does
            # not actually exist. Just pass, removing the contact
            # from the current contact set.
            except:
                parContacts.remove(contact);
    return edgeSet;