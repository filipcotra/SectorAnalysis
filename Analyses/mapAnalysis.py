import json;

# Defining constants for edgeSet value indices.
i_PAR_NAME = 0;
i_PAR_RES = 1;
i_CONTACT_NAME = 3;
i_CONTACT_RES = 4;

# Purpose: To create contact maps from a given set of edges.
# Parameters:
#   edgeSet = The set of edges describing the structure.
#   aaSeq = The amino-acid sequence of the structure.
#   pdb = The pdb code for the structure.
#   outputDir = The directory where output will be stored.
def createMaps(edgeSet, aaSeq, pdb, outputDir):
    seqList = [];
    for i in aaSeq:
        seqList.append(aaSeq[i]);
    seqStr = "".join(seqList);
    outputName = f"{outputDir}/{pdb}_CM.txt";
    # Making headers.
    outputLine = f"##\t{seqStr}\n#\tExperimental Map - {pdb}\n";
    # Storing maps in dictionaries.
    BB2BB_map = {};
    SC2SC_map = {};
    SC2BB_map = {};
    # Iterating through the edges of the edge set, populating the dictionaries.
    for edge in edgeSet:
        # Extracting info from the sets.
        srcPar = edge[i_PAR_NAME];
        srcRes = edge[i_PAR_RES];
        dstPar = edge[i_CONTACT_NAME];
        dstRes = edge[i_CONTACT_RES];
        # Finding AA types.
        srcType = srcPar[-1];
        dstType = dstPar[-1];
        # Populating contact map dictionaries.
        if srcType == "0":
            if dstType == "0": # BB2BB. Adding both directions.
                addVal(srcRes, dstRes, BB2BB_map); # Adding dstRes to srcRes key.
                addVal(dstRes, srcRes, BB2BB_map); # Adding srcRes to dstRec key.
            else: # SC2BB. Adding srcRes to dstRes.
                addVal(dstRes, srcRes, SC2BB_map);
        else: # srcType == 1.
            if dstType == "0": # SC2BB. Adding dstRes to srcRes.
                addVal(srcRes, dstRes, SC2BB_map);
            else: # SC2SC. Adding both directions.
                addVal(srcRes, dstRes, SC2SC_map); # Adding dstRes to srcRes.
                addVal(dstRes, srcRes, SC2SC_map); # Adding srcRes to dstRes.
    # Outputting in json format.
    outputLine += f"{json.dumps(BB2BB_map)}\n{json.dumps(SC2SC_map)}\n{json.dumps(SC2BB_map)}\n"
    with open(outputName, "w") as outFile:
        outFile.write(outputLine);

# Purpose: To add a value to a dictionary.
# Parameters:
#   key = The key value.
#   val = The value.
#   df = The dataframe.
def addVal(key, val, df):
    if key in df:
        if val not in df[key]:
            df[key].append(val);
    else:
        df[key] = [val];

