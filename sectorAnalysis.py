import sys;
from Analyses.pairwiseAnalysis import analyze, export;
from Core.fetchInfo import fetchLineInfo, filterSet;
from Core.extractEdgeSet import getEdgeSet;

# Defining constants for fetchInfo returns.
i_PAR_NAME = 0;
i_PAR_RES = 3;
i_EXTRA_CONTACTS = 5;
# Only input should be a file with a list of files.
files = open(sys.argv[sys.argv.index("-f") + 1], "r");
outputSuff = sys.argv[sys.argv.index("-o") + 1];
outputDir = sys.argv[sys.argv.index("-d") + 1];
# Tracking integer IDs representing distinct PDB stuctures and
# particles.
structID = 0;
parID = 0;
SOLV_PAR = ("O0", "-");

# Looping through the files and tracking the contact set of each particle.
for fileName in files:
    parIDMap = {SOLV_PAR: parID}; # Initializing the solvent molecule as residue 0.
    parID += 1;
    # Getting information from the file.
    print(f"Analyzing: {fileName}");
    currFile = open(fileName.rstrip(), "r");
    # Tracking the contact set of each particle for the current file.
    particleContacts = {SOLV_PAR: []};
    # In this file, identifying the sectors for non-backbone contacts and
    # adding the respective amino acids to the sector values. Fetching info
    # by particle.
    lastLine = False;
    currPar = None;
    lastPar = None;
    contactSet = []; # Want to track all the contact info for the current particle.
    while not lastLine:
        line = currFile.readline();
        # Checking if we are at the end of the line.
        currentPos = currFile.tell(); # Identifying the current position of the file.
        lastLine = True if currFile.readline() == "" else False;
        currFile.seek(currentPos); # Resetting the file position to the current line.
        # If we run into any exceptions, which should be rare, just
        # continue to the next.
        try:
            lineInfo = fetchLineInfo(line);
        except:
            lineInfo = False;
        # If lineInfo is False, continue to the next line.
        # However, if we have reached the last line, don't.
        if not lineInfo and not lastLine:
            continue;
        # As long as lineInfo isn't false, we should collect new information.
        if lineInfo != False:
            # Importantly, this key uses the int version of the
            # particle residue number.
            currPar = (lineInfo[i_PAR_NAME], lineInfo[i_PAR_RES]);
            if lastPar is None:
                lastPar = currPar;
            # If the current particle is different from the last, push the
            # current contact set.
            if currPar != lastPar:
                filteredSet = filterSet(contactSet);
                # Noting the solvent contacts.
                solvContacts = [x for x in filteredSet if "O0" in x];
                revSolv = [x[4:8] + x[0:4] + (x[8],) + (-1,) + (x[10],) for x in solvContacts];
                particleContacts[SOLV_PAR] += revSolv;
                # Noting the contact set and then resetting it to empty.
                particleContacts[lastPar] = filteredSet;
                parIDMap[lastPar] = parID;
                parID += 1;
                contactSet = [];
                lastPar = currPar;
            # Always adding the current line to the new contact set.
            contactSet.append(lineInfo);
        # If we are on the last line, push the current contact set.
        if lastLine:
            filteredSet = filterSet(contactSet);
            solvContacts = [x for x in filteredSet if "O0" in x];
            revSolv = [x[4:8] + x[0:4] + (x[8],) + (-1,) + (x[10],) for x in solvContacts];
            particleContacts[SOLV_PAR] += revSolv;
            particleContacts[currPar] = filteredSet;
            parIDMap[currPar] = parID;
            parID += 1;
    currFile.close();
    # Now we have all the contact sets for each particle
    # in the file, which we should extract to a set of overall
    # edges for the dataset.
    edgeSet = getEdgeSet(particleContacts);
    # Pairwise analysis will look at pairwise probabilities, include sector-sector
    # contacts, contacts between amino-acid types, and the number of shared contacts.
    # This will require both the total edge set and the contact set for each particle.
    analyze(edgeSet, particleContacts, parIDMap, structID);
    # Deleting some of the bigger objects to hopefully save memory.
    del contactSet, particleContacts, edgeSet;
    structID += 1; # Iterating the structure ID each time a new structure file is opened.
files.close();
export(outputDir, outputSuff);
