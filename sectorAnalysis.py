import sys;
from Core.fetchInfo import fetchLineInfo;
from Analyses.particleCapacity import capacityAnalysis, printParticleOutput;
from Analyses.sectorCount import countAnalysis, printCountOutput;
from Analyses.sectorOccupancy import occupancyAnalysis, printOccupancyOutput;

# Only input should be a file with a list of files.
files = open(sys.argv[sys.argv.index("-f") + 1], "r");
# Checking which analyses will be performed.
performParticle = int(sys.argv[sys.argv.index("-p") + 1]);
performCount = int(sys.argv[sys.argv.index("-c") + 1]);
performOccupancy = int(sys.argv[sys.argv.index("-o") + 1]);

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

# Looping through the files, calling the necessary analysis modules
# as necessary, per line.
for fileName in files:
    print(f"Analyzing: {fileName}");
    currFile = open(fileName.rstrip(), "r");
    # In this file, identifying the sectors for non-backbone contacts and
    # adding the respective amino acids to the sector values.
    for line in currFile:
        lineInfo = fetchLineInfo(line);
        # If fetchLineInfo returns False, continue to the next line.
        if not lineInfo:
            continue;
        # Call appropriate functions.
        if performParticle:
            capacityAnalysis(lineInfo[i_PAR_NAME], lineInfo[i_PAR_RES]);
        if performCount:
            countAnalysis(lineInfo[i_PAR_NAME], lineInfo[i_PAR_TYPE], lineInfo[i_CONTACT_NAME],
                          lineInfo[i_CONTACT_TYPE], lineInfo[i_RES_DIFF], lineInfo[i_SECTOR_NUM]);
        if performOccupancy:
            occupancyAnalysis(lineInfo[i_PAR_NAME], lineInfo[i_PAR_RES], lineInfo[i_CONTACT_NAME],
                              lineInfo[i_CONTACT_RES], lineInfo[i_SECTOR_NUM]);
    currFile.close();
# Printing all output.
printParticleOutput();
printCountOutput();
printOccupancyOutput();

