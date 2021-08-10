
import os
import re
import sys
import suffer

digPat = re.compile("\\d+")

protLoc = sys.argv[1]
workLoc = sys.argv[2]
lookFile = sys.argv[3]
profLoc = sys.argv[4]
lookTSV = sys.argv[5]

suffArr = suffer.ProteosSuffixWrapper(profLoc, protLoc)
pepSearch = suffer.ProteosPeptideSearch(lookFile, workLoc)

# load in the peptide sequences
lookTSVF = open(lookTSV, "r")
tsvFL = True
pepSeqs = []
for line in lookTSVF:
    if line.strip() == "":
        continue
    if tsvFL:
        tsvFL = False
        continue
    lineS = line.split()
    pepSeqs.append(lineS[0])
lookTSVF.close()
curPSI = 0

# load ENSP/G/T data from protLoc
varIDF = open(os.path.join(protLoc, "varids.tsv"))
pDatMap = {}
for line in varIDF:
    if line.strip() == "":
        continue
    lineS = line.split()
    protNam = lineS[0]
    perInd = protNam.find(".")
    if perInd >= 0:
        protNam = protNam[:perInd]
    tranNam = lineS[1]
    perInd = tranNam.find(".")
    if perInd >= 0:
        tranNam = tranNam[:perInd]
    geneNam = lineS[2]
    perInd = geneNam.find(".")
    if perInd >= 0:
        geneNam = geneNam[:perInd]
    pDatMap[protNam] = [tranNam, geneNam]
varIDF.close()

outLoc = os.path.join(workLoc, "proteos.tsv")
outF = open(outLoc, "w")
toPrin = []
toPrin.append("peptide_seq") #The sequence of the peptide
toPrin.append("peptide_seq_ll") #The sequence of the peptide with isoleucine converted to leucine.
toPrin.append("peptide_hit") #Whether this current entry is a hit.
toPrin.append("ensembl_protein_id") #The Ensemble ID of a protein the peptide can be found in.
toPrin.append("ensembl_transcript_id") #The corresponding transcript ID.
toPrin.append("ensemble_gene_id") #The corresponding gene ID.
toPrin.append("chromosome") #The chromosome the protein is from.
toPrin.append("protein_location_genome") #The location of the originating bases of the variant of the protein the peptide can be found in.
toPrin.append("peptide_location_genome") #The bases in the genome the peptide originates from.
toPrin.append("peptide_start_transcript") #The first amino acid the peptide corresponds to in the variant.
toPrin.append("peptide_end_transcript") #The last amino acid the peptide corresponds to in the variant.
toPrin.append("peptide_start_reference") #The first amino acid in the reference corresponding to a part of the peptide.
toPrin.append("peptide_end_reference") #The last amino acid in the reference corresponding to a part of the peptide.
toPrin.append("snp_location_genome") #The nucleotide changes that lead to the current variant of the protein the peptide matches.
toPrin.append("sap_location_reference") #The corresponding sap changes.
toPrin.append("sap_location_transcript") #The corresponding locations in the variant protein.
toPrin.append("sap_location_peptide") #The locations of any saps within the peptide.
outF.write("\t".join(toPrin) + "\n")
numCols = len(toPrin)

curEnt = pepSearch.getNextMatch()
while not (curEnt is None):
    # need to get
    # peptide_seq_
    # peptide_seq_ll
    # peptide_hit
    peptideSeqL = curEnt.sequence
    peptideSeq = pepSeqs[curPSI]
    curPSI = curPSI + 1
    if len(curEnt.matches) == 0:
        outF.write(peptideSeq + "\t" + peptideSeqL + "\t" + "FALSE" + ("\t"*(numCols-3)) + "\n")
        curEnt = pepSearch.getNextMatch()
        continue
    for mi in range(len(curEnt.names)):
        matDat = curEnt.matches[mi]
        tranCig = curEnt.varCigs[mi]
        locData = curEnt.protLocs[mi]
        # ensembl_protein_id
        # ensembl_transcript_id
        # ensemble_gene_id
        ensGeneID = ""
        ensTranID = ""
        ensProtID = curEnt.names[mi]
        undInd = ensProtID.find("_")
        if undInd >= 0:
            ensProtID = ensProtID[:undInd]
        if ensProtID in pDatMap:
            curPDat = pDatMap[ensProtID]
            ensTranID = curPDat[0]
            ensGeneID = curPDat[1]
        # snp_location_genome
        # sap_location_reference
        snpLocGen = ""
        sapLocPro = ""
        sapLocations = []
        aaChange = curEnt.varDiff[mi]
        for aaCSet in aaChange:
            for sinChange in aaCSet.split(';'):
                aanaSpl = sinChange.split(':')
                if len(aanaSpl) != 2:
                    sys.stderr.write(ensProtID + ": Questionable variant " + sinChange + "\n")
                digM = digPat.search(aanaSpl[0])
                if not (digM is None):
                    sapLocations.append(int(aanaSpl[0][digM.span()[0] : digM.span()[1]]))
                sapLocPro = sapLocPro + ("," if (len(sapLocPro) > 0) else "") + aanaSpl[0]
                snpLocGen = snpLocGen + ("," if (len(snpLocGen) > 0) else "") + aanaSpl[1]
        sapLocations = set(sapLocations)
        # chromosome
        # protein_location_genome
        # peptide_location_genome
        # peptide_start_transcript
        # peptide_end_transcript
        # peptide_start_reference
        # peptide_end_reference
        matStart = matDat[2]
        matEnd = matDat[3]
        pepStartTran = repr(1 + matStart)
        pepEndTran = repr(matEnd)
        proLocGen = ""
        pepLocGen = ""
        pepStartPro = ""
        pepEndPro = ""
        chrom = ""
        sapLocTran = ""
        sapLocPep = ""
        if len(locData) != 0:
            chrom = locData[0][0]
            # expand the location data into triples
            isPos = locData[0][3]
            haveWarn = False
            fullLocs = []
            for ldat in locData:
                if ldat[3] != isPos:
                    sys.stderr.write(ensProtID + ": Transcript uses both strands...\n")
                    break
                for i in range(ldat[1], ldat[2]):
                    fullLocs.append(i)
                    if (ldat[0] != chrom) and not haveWarn:
                        sys.stderr.write(ensProtID + ": Multiple chromosomes in location data " + chrom + " and " + ldat[0] + "\n")
            if not isPos:
                fullLocs = list(reversed(fullLocs))
            tmpLocs = []
            i = 0
            while i < len(fullLocs):
                if i+3 > len(fullLocs):
                    sys.stderr.write(ensProtID + ": location data not divisible by 3.\n")
                    break
                tmpLocs.append([fullLocs[i], fullLocs[i+1], fullLocs[i+2]])
                i = i + 3
            fullLocs = tmpLocs
            # figure out the mapping for the cigars
            tmpCig = [i for i in range(matEnd)]
            if tranCig == "???":
                sys.stderr.write(ensProtID + ": missing cigar")
            else:
                try:
                    cigReg = suffer.cigarStringToReferencePositions(0, tranCig)
                    if (cigReg[0] > 0) or (cigReg[1] > 0):
                        raise ValueError("Soft clip in protein cigar.")
                    if len(cigReg[2]) < matEnd:
                        raise ValueError("Search reports match beyond end of the cigar.")
                    tmpCig = cigReg[2]
                except ValueError as err:
                    errMess = getattr(err, 'message', repr(err))
                    sys.stderr.write(ensProtID + ": bad cigar " + tranCig + ": " + errMess + "\n")
            # make sure fullLocs has the right length
            haveWarnCigGff = False
            for i in range(len(tmpCig)):
                curCI = tmpCig[i]
                if curCI >= len(fullLocs):
                    if not haveWarnCigGff:
                        haveWarnCigGff = True
                        sys.stderr.write(ensProtID + ": cigar extends beyond location data\n")
                    tmpCig[i] = -1
            # protein_location_genome: chr:start:end:direction  1based location of the start of the protein
            proLocGen = []
            for ldat in locData:
                if len(proLocGen) > 0:
                    proLocGen.append("|")
                proLocGen.extend(":".join([ldat[0], repr(ldat[1]), repr(ldat[2]), "+" if ldat[3] else "-"]))
            proLocGen = "".join(proLocGen)
            # peptide_location_genome: chr:start-end  same as before
            pepLocs = []
            for i in range(matStart,matEnd):
                curRI = tmpCig[i]
                if curRI < 0:
                    pepLocs.append(-1)
                else:
                    pepLocs.extend(fullLocs[curRI])
            expDiff = 1 if isPos else -1
            i = 0
            while i < len(pepLocs):
                if pepLocs[i] < 0:
                    while i < len(pepLocs):
                        if pepLocs[i] >= 0:
                            break
                        i = i + 1
                    pepLocGen = pepLocGen + ("|" if (len(pepLocGen) > 0) else "") + "-1"
                else:
                    pepLocGen = pepLocGen + ("|" if (len(pepLocGen) > 0) else "") + chrom + ":" + repr(pepLocs[i])
                    origLoc = pepLocs[i]
                    lastLoc = pepLocs[i]
                    i = i + 1
                    while i < len(pepLocs):
                        if pepLocs[i] != (lastLoc + expDiff):
                            break
                        lastLoc = pepLocs[i]
                        i = i + 1
                    if origLoc != lastLoc:
                        pepLocGen = pepLocGen + "-" + repr(lastLoc)
            # peptide_start_reference: 1 based index of the start of the peptide in the reference protein
            # peptide_end_reference: 1 based index of the end of the peptide in the reference protein
            pepStartPro = len(tmpCig)
            pepEndPro = -1
            for i in range(matStart,matEnd):
                cigI = tmpCig[i]
                if cigI >= 0:
                    pepStartPro = min(pepStartPro, cigI)
                    pepEndPro = max(pepEndPro, cigI)
            pepStartPro = repr(pepStartPro+1)
            pepEndPro = repr(pepEndPro+1)
            # find any locations in the transcript that are also in sapLocations
            for i in range(len(tmpCig)):
                if tmpCig[i] in sapLocations:
                    sapLocTran = sapLocTran + (";" if (len(sapLocTran) > 0) else "") + repr(i)
                    if (i >= matStart) and (i < matEnd):
                        sapLocPep = sapLocPep + (";" if (len(sapLocPep) > 0) else "") + repr(i - matStart)
        # print the results
        toPrin = []
        toPrin.append(peptideSeq)
        toPrin.append(peptideSeqL)
        toPrin.append("TRUE")
        toPrin.append(ensProtID)
        toPrin.append(ensTranID)
        toPrin.append(ensGeneID)
        toPrin.append(chrom)
        toPrin.append(proLocGen)
        toPrin.append(pepLocGen)
        toPrin.append(pepStartTran)
        toPrin.append(pepEndTran)
        toPrin.append(pepStartPro)
        toPrin.append(pepEndPro)
        toPrin.append(snpLocGen)
        toPrin.append(sapLocPro)
        toPrin.append(sapLocTran)
        toPrin.append(sapLocPep)
        outF.write("\t".join(toPrin) + "\n")
    curEnt = pepSearch.getNextMatch()


