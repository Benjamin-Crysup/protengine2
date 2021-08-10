
import os
import sys
import suffer

protLoc = sys.argv[1]
workLoc = sys.argv[2]
lookFile = sys.argv[3]
profLoc = sys.argv[4]

suffArr = suffer.ProteosSuffixWrapper(profLoc, protLoc)

for i in range(5, len(sys.argv)):
    phase = int(sys.argv[i])
    if phase == 1:
        suffer.runPeptideSearchSearch(suffArr, lookFile, workLoc, sys.stderr)
    elif phase == 2:
        suffer.runPeptideSearchNames(suffArr, lookFile, workLoc, sys.stderr)
    elif phase == 3:
        suffer.runPeptideSearchProtein(suffArr, lookFile, workLoc, sys.stderr)
    elif phase == 4:
        suffer.runPeptideSearchVariant(suffArr, lookFile, workLoc, sys.stderr)
    elif phase == 0:
        suffer.runPeptideSearch(suffArr, lookFile, workLoc, sys.stderr)
    elif phase == 5:
        outLoc = os.path.join(workLoc, "freqs.tsv")
        pepSearch = suffer.ProteosPeptideSearch(lookFile, workLoc)
        outF = open(outLoc, "w")
        curEnt = pepSearch.getNextMatch()
        while not (curEnt is None):
            outF.write(curEnt.sequence + "\t" + repr(len(curEnt.getFullPeptideOwnership())) + "\t" + "|".join(curEnt.getFullPeptideVariants()) + "\n")
            curEnt = pepSearch.getNextMatch()
        outF.close()

