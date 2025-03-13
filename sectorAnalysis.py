import sys;
from Core.fetchInfo import fetchLineInfo, filterSet;
from Analyses.particleCapacity import capacityAnalysis, printParticleOutput;
from Analyses.sectorCount import countAnalysis, printCountOutput;
from Analyses.sectorOccupancy import occupancyAnalysis, printOccupancyOutput;

# Defining constants for fetchInfo returns.
i_PAR_NAME = 0;
i_PAR_RES = 3;
# Only input should be a file with a list of files.
files = open(sys.argv[sys.argv.index("-f") + 1], "r");
# Checking which analyses will be performed.
performCapacity = int(sys.argv[sys.argv.index("-p") + 1]);
performCount = int(sys.argv[sys.argv.index("-c") + 1]);
performOccupancy = int(sys.argv[sys.argv.index("-o") + 1]);

# Looping through the files, calling the necessary analysis modules
# as necessary, per line.
for fileName in files:
    print(f"Analyzing: {fileName}");
    currFile = open(fileName.rstrip(), "r");
    # In this file, identifying the sectors for non-backbone contacts and
    # adding the respective amino acids to the sector values. Fetching info
    # by particle.
    currPar = None;
    lastPar = None;
    contactSet = []; # Want to track all the contact info for the current particle.
    for line in currFile:
        lineInfo = fetchLineInfo(line);
        # If fetchLineInfo returns False, continue to the next line.
        if not lineInfo:
            continue;
        currPar = (lineInfo[i_PAR_NAME], lineInfo[i_PAR_RES]);
        if lastPar is None:
            lastPar = currPar;
        # If we encounter a new particle, call the analysis functions.
        if lastPar != currPar:
            filteredSet = filterSet(contactSet);
            # Call appropriate functions. contactSet will only contain
            # real contacts.
            if performCapacity:
                capacityAnalysis(lastPar, filteredSet);
            if performCount:
                countAnalysis(filteredSet);
            if performOccupancy:
                occupancyAnalysis(filteredSet);
            contactSet = [];
            lastPar = currPar;
        # Adding the contact to the contact set regardless.
        contactSet.append(lineInfo);
    currFile.close();
# Printing all output.
printParticleOutput();
printCountOutput();
printOccupancyOutput();

