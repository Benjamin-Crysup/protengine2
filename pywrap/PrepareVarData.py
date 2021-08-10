
import sys

# Prepare a TSV containing data on variants.
# Needs a lookup table containing which individuals have which variant: assumed to be sorted on variant name

fromLUT = open(sys.argv[1], "r")
lutLin = fromLUT.readline()
lutLin = fromLUT.readline()

while any([lutLin != ""]):
    curProN = None
    cigPosibs = {}
    diffPosibs = set()
    curVInds = []
    # get data from the lookup table
    while lutLin != "":
        if lutLin.strip() == "":
            lutLin = fromLUT.readline()
            continue
        lutS = lutLin.split()
        if (not (curProN is None)) and (lutS[0] != curProN):
            break
        curProN = lutS[0]
        if lutS[2] in cigPosibs:
            cigPosibs[lutS[2]] = cigPosibs[lutS[2]] + 1
        else:
            cigPosibs[lutS[2]] = 1
        if len(lutS) > 3:
            diffPosibs.add(lutS[3])
        curVInds.append(lutS[1])
        #TODO the vcf data
        lutLin = fromLUT.readline()
    # dump out for the variant
    if curProN is None:
        continue
    if len(cigPosibs) == 0:
        raise IOError("Missing CIGAR for " + curProN)
    if len(cigPosibs) > 1:
        sys.stderr.write("Inconsitent CIGAR for " + curProN + "\n");
    bestC = 0
    curVCig = None
    for cig in cigPosibs:
        if cigPosibs[cig] > bestC:
            bestC = cigPosibs[cig]
            curVCig = cig
    toPrin = [curProN, curVCig, "INDS{"] + curVInds + ["}SDNI", "DIFS{"] + list(diffPosibs) + ["}SFID"]
    print("\t".join(toPrin))

fromLUT.close()
