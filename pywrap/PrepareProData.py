
import sys

# Prepare a TSV containing data on proteins.
# Needs a bed file containing location information: assumed to be sorted on protein name

fromBed = open(sys.argv[1], "r")
bedLin = fromBed.readline()

while any([bedLin != ""]):
    curProN = None
    curProLocD = []
    # get data from the lookup table
    while bedLin != "":
        if bedLin.strip() == "":
            bedLin = fromBed.readline()
            continue
        lutS = bedLin.split()
        if (not (curProN is None)) and (lutS[3] != curProN):
            break
        curProN = lutS[3]
        curProLocD.append(lutS[0])
        curProLocD.append(lutS[1])
        curProLocD.append(lutS[2])
        curProLocD.append(lutS[4])
        bedLin = fromBed.readline()
    # dump out for the variant
    if curProN is None:
        continue
    toPrin = [curProN, "LOCS{"] + curProLocD + ["}SCOL"]
    print("\t".join(toPrin))

fromBed.close()



