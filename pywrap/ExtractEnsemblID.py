'''Extracts Ensembl IDs from a fasta file.'''

import re
import sys

patP = re.compile("ENSP\\d+\\.\\d+")
patT = re.compile("ENST\\d+\\.\\d+")
patG = re.compile("ENSG\\d+\\.\\d+")

if (len(sys.argv) == 1) or (sys.argv[1] == "-"):
    fromF = sys.stdin
else:
    fromF = open(sys.argv[1])

print("ProteinID\tTranscriptID\tGeneID")
for line in fromF:
    if line.strip() == "":
        continue
    if line.strip()[0] != ">":
        continue
    matP = patP.search(line)
    matT = patT.search(line)
    matG = patG.search(line)
    if matP is None:
        matP = ""
    else:
        matP = line[matP.start():matP.end()]
    if matT is None:
        matT = ""
    else:
        matT = line[matT.start():matT.end()]
    if matG is None:
        matG = ""
    else:
        matG = line[matG.start():matG.end()]
    print(matP + "\t" + matT + "\t" + matG)
fromF.close()

