import os
import re
import sys
import struct
import subprocess

#reference needed files
#seqs.gail          : packaged sequence data
#protdata.tsv       : data on the proteins
#vardata.tsv        : data on the variants
#digest.dig         : digest specification
#people.tsv         : all the people/haplotypes in the database
#
#optional files
#seqs.gail.sa       : suffix array
#protdata.bctsv     : data on the proteins, sorted on first column (protein name)
#vardata.bctsv      : data on the variants, sorted on first column (variant name)
#varids.tsv         : all the protein variants in the database
#

def cigarStringToReferencePositions(refPos0, cigStr):
    fillPos = []
    curRef = refPos0
    seenAction = False
    softClipStart = 0
    softClipEnd = 0
    seenHardEnd = False
    ci = 0
    while ci < len(cigStr):
        cd = ci
        while cd < len(cigStr):
            if not (cigStr[cd] in "+-0123456789"):
                break
            cd = cd + 1
        if cd == ci:
            raise ValueError("Cigar operation missing count.")
        if cd == len(cigStr):
            raise ValueError("Cigar operation missing operation.")
        opCount = int(cigStr[ci:cd])
        if opCount < 0:
            raise ValueError("Negative operation count.")
        opChar = cigStr[cd]
        if (opChar == "M") or (opChar == "=") or (opChar == "X"):
            if (softClipEnd > 0) or seenHardEnd:
                raise ValueError("Match operation in cigar string after clipping.")
            seenAction = True
            for i in range(opCount):
                fillPos.append(curRef)
                curRef = curRef + 1
        elif (opChar == "I"):
            if (softClipEnd > 0) or seenHardEnd:
                raise ValueError("Insertion operation in cigar string after clipping.")
            seenAction = True
            for i in range(opCount):
                fillPos.append(-1)
        elif (opChar == "D") or (opChar == "N"):
            if (softClipEnd > 0) or seenHardEnd:
                raise ValueError("Deletion operation in cigar string after clipping.")
            seenAction = True
            curRef = curRef + opCount
        elif (opChar == "S"):
            if seenHardEnd:
                raise ValueError("Cannot soft clip after hard clip.")
            if seenAction:
                softClipEnd = softClipEnd + opCount
            else:
                softClipStart = softClipStart + opCount
        elif (opChar == "H"):
            if seenAction or (softClipStart > 0):
                seenHardEnd = True
        elif (opChar == "P"):
            # do not care
            pass
        else:
            raise ValueError("Unknown cigar operation.")
        ci = cd + 1
    return [softClipStart, softClipEnd, fillPos]

class ProteosSuffixWrapper:
    """A suffix array to search through."""
    def __init__(self, progPath, protFolder):
        """
            Wrap a suffix array.
            @param progFolder: The path to the profinman program.
            @param protFolder: The folder containing the protein database.
        """
        self.progPath = progPath
        """The path to the profinman program."""
        self.protFolder = protFolder
        """The folder containing the protein database."""
    def listAllPeople(self):
        '''
            Get all the people in this database.
            @return: The set of people.
        '''
        pplFN = os.path.join(self.protFolder, "people.tsv")
        if not os.path.exists(pplFN):
            raise IOError("Missing " + pplFN)
        allPPL = []
        pplF = open(pplFN)
        for line in pplF:
            if line.strip() == "":
                continue
            allPPL.append(line.strip())
        pplF.close()
        return allPPL
    #TODO get info on the thing


def runPeptideSearchSearch(wrapped, lookFor, workFold, repStr):
    # basic setup
    protFolder = wrapped.protFolder
    progPath = wrapped.progPath
    refFileN = os.path.join(protFolder, "seqs.gail")
    refSufAN = os.path.join(protFolder, "seqs.gail.sa")
    os.makedirs(os.path.abspath(workFold), exist_ok=True)
    # set up the stream: make a temporary file
    lookForFa = None
    searchStr = None
    if isinstance(lookFor, str):
        lookForFa = lookFor
    else:
        lookForFa = os.path.join(workFold, "searchList.fa")
        tmpFAO = open(lookForFa, "w")
        for pep in lookFor:
            tmpFAO.write(">lookie\n" + pep + "\n")
        tmpFAO.close()
    searchStr = open(lookForFa, "rb")
    # perform the search (and sort the results)
    repStr.write("Running search...")
    repStr.flush()
    searchResFN = os.path.join(workFold, "searchInds.bin")
    searchProc = None
    if os.path.exists(refSufAN):
        searchProc = subprocess.Popen([progPath, "findsa", "--ref", refFileN, "--sa", refSufAN], stdin=searchStr, stdout=subprocess.PIPE)
    else:
        searchProc = subprocess.Popen([progPath, "findfa", "--ref", refFileN], stdin=searchStr, stdout=subprocess.PIPE)
    stripProc = subprocess.Popen([progPath, "matfildig", "--dig", os.path.join(protFolder, "digest.dig"), "--ref", refFileN], stdin=searchProc.stdout, stdout=subprocess.PIPE)
    sortfProc = subprocess.Popen([progPath, "sortfound", "--work", workFold, "--out", searchResFN], stdin=stripProc.stdout)
    if searchProc.wait() != 0:
        raise IOError("Problem running the initial search.")
    if stripProc.wait() != 0:
        raise IOError("Problem filtering inconsistent digests.")
    if sortfProc.wait() != 0:
        raise IOError("Problem sorting the initial search.")
    searchStr.close()
    repStr.write("DONE\n")

def runPeptideSearchNames(wrapped, lookFor, workFold, repStr):
    progPath = wrapped.progPath
    protFolder = wrapped.protFolder
    refFileN = os.path.join(protFolder, "seqs.gail")
    searchResFN = os.path.join(workFold, "searchInds.bin")
    # get the names of the proteins
    repStr.write("Getting names...")
    repStr.flush()
    nameResFN = os.path.join(workFold, "searchNames.tsv")
    getnProc = subprocess.Popen([progPath, "nameget", "--out", nameResFN, "--ref", refFileN, "--match", searchResFN])
    if getnProc.wait() != 0:
        raise IOError("Problem getting match names.")
    repStr.write("DONE\n")

def loadInTSV(filName):
    tFil = open(filName, "r")
    allLines = {}
    for line in tFil:
        if line.strip() == "":
            continue
        lineS = line.split()
        if lineS[0] in allLines:
            raise IOError("Duplicate entries in database: " + lineS[0])
        allLines[lineS[0]] = line
    tFil.close()
    return allLines

def runPeptideSearchProtein(wrapped, lookFor, workFold, repStr):
    progPath = wrapped.progPath
    protFolder = wrapped.protFolder
    nameResFN = os.path.join(workFold, "searchNames.tsv")
    # lookup the protein data (and sort)
    repStr.write("Getting protein data...")
    repStr.flush()
    protDataFN = os.path.join(workFold, "protdata.tsv")
    if os.path.exists(os.path.join(protFolder, "protdata.bctsv")):
        lookPProc = subprocess.Popen([progPath, "findtabs", "--dcol", "0", "--fcol", "0", "--in", os.path.join(protFolder, "protdata.bctsv")], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        sortPProc = subprocess.Popen([progPath, "sortfindtab", "--dump", protDataFN, "--work", os.path.normpath(workFold)], stdin=lookPProc.stdout)
        nameResF = open(nameResFN, "r")
        for line in nameResF:
            if line.strip() == "":
                continue
            lineS = line.strip().split("_")
            lookPProc.stdin.write((lineS[0] + "\n").encode())
        nameResF.close()
        lookPProc.stdin.close()
        if lookPProc.wait() != 0:
            raise IOError("Problem searching for protein data.")
        if sortPProc.wait() != 0:
            raise IOError("Problem sorting protein data.")
    elif os.path.exists(os.path.join(protFolder, "protdata.tsv")):
        fullDBase = loadInTSV(os.path.join(protFolder, "protdata.tsv"))
        curInd = 0
        outF = open(protDataFN, "w")
        nameResF = open(nameResFN, "r")
        for line in nameResF:
            if line.strip() == "":
                continue
            lineS = line.strip().split("_")
            if lineS[0] in fullDBase:
                outF.write(repr(curInd) + "\t" + fullDBase[lineS[0]])
            curInd = curInd + 1
        outF.close()
        nameResF.close()
    else:
        raise IOError("No protein database found: need protdata.tsv or protdata.bctsv.")
    repStr.write("DONE\n")

def runPeptideSearchVariant(wrapped, lookFor, workFold, repStr):
    progPath = wrapped.progPath
    protFolder = wrapped.protFolder
    nameResFN = os.path.join(workFold, "searchNames.tsv")
    # lookup the variant data (and sort)
    repStr.write("Getting variant data...")
    repStr.flush()
    varDataFN = os.path.join(workFold, "vardata.tsv")
    if os.path.exists(os.path.join(protFolder, "vardata.bctsv")):
        nameResF = open(nameResFN, "rb")
        lookVProc = subprocess.Popen([progPath, "findtabs", "--dcol", "0", "--fcol", "0", "--in", os.path.join(protFolder, "vardata.bctsv")], stdout=subprocess.PIPE, stdin=nameResF)
        sortVProc = subprocess.Popen([progPath, "sortfindtab", "--dump", varDataFN, "--work", os.path.normpath(workFold)], stdin=lookVProc.stdout)
        if lookVProc.wait() != 0:
            raise IOError("Problem searching for protein variant data.")
        if sortVProc.wait() != 0:
            raise IOError("Problem sorting protein variant data.")
        nameResF.close()
    elif os.path.exists(os.path.join(protFolder, "vardata.tsv")):
        fullDBase = loadInTSV(os.path.join(protFolder, "vardata.tsv"))
        curInd = 0
        outF = open(varDataFN, "w")
        nameResF = open(nameResFN, "r")
        for line in nameResF:
            if line.strip() == "":
                continue
            lineS = line.strip()
            if lineS in fullDBase:
                outF.write(repr(curInd) + "\t" + fullDBase[lineS])
            curInd = curInd + 1
        outF.close()
        nameResF.close()
    else:
        raise IOError("No variant database found: need vardata.tsv or vardata.bctsv.")
    repStr.write("DONE\n")

def runPeptideSearch(wrapped, lookFor, workFold, repStr):
    runPeptideSearchSearch(wrapped, lookFor, workFold, repStr)
    runPeptideSearchNames(wrapped, lookFor, workFold, repStr)
    runPeptideSearchProtein(wrapped, lookFor, workFold, repStr)
    runPeptideSearchVariant(wrapped, lookFor, workFold, repStr)

_digPat = re.compile("\\d+")

class ProteosParsedDifference:
    """A difference."""
    # 167-168S>F:109130043C>T
    def __init__(self, difStr):
        """
            Parse a difference.
        """
        self.warning = None
        """A warning while parsing this difference."""
        self.protRange = None
        """The range the difference covers in the reference."""
        self.protChange = ""
        """A string describing the protein change."""
        self.genLoc = 0
        """The genome location of the change."""
        self.genChange = ""
        """A string describing the genetic change."""
        self.diffDesc = difStr
        """The full difference string"""
        aanaSpl = difStr.split(':')
        if len(aanaSpl) != 2:
            self.warning = "Questionable variant: " + difStr
            return
        protDs = aanaSpl[0].split("-")
        if len(protDs) != 2:
            self.warning = "Malformed amino acid variant data: " + difStr
            return
        protSs = protDs[0]
        protEs = protDs[1]
        genD = aanaSpl[1]
        digMss = _digPat.search(protSs)
        if digMss is None:
            self.warning = "Missing amino acid change start."
            return
        digMes = _digPat.search(protEs)
        if digMes is None:
            self.warning = "Missing amino acid change end."
            return
        self.protRange = [0,0]
        self.protRange[0] = int(protSs[digMss.span()[0] : digMss.span()[1]])
        self.protRange[1] = int(protEs[digMes.span()[0] : digMes.span()[1]])
        self.protChange = protEs[digMes.span()[1]:]
        digMps = _digPat.search(genD)
        if digMps is None:
            self.warning = "Missing genetic change location."
            return
        self.genLoc = int(genD[digMps.span()[0] : digMps.span()[1]])
        self.genChange = genD[digMps.span()[1]:]
    def __hash__(self):
        return hash(self.diffDesc)
    def __eq__(self, other):
        return hasattr(other, "diffDesc") and (self.diffDesc == other.diffDesc)
    def __repr__(self):
        return self.diffDesc

def parseProteosDifferences(fullDiffStr):
    """
        Parse all the differences relevant to a protein variant.
        @param fullDiffStr: The full difference string.
        @return: A (possibly empty) list of differences.
    """
    return [ProteosParsedDifference(sinC) for sinC in fullDiffStr.split(";")]

class ProteosPeptideInfo:
    """Search result info for a peptide."""
    def __init__(self, seq, mat, nam, pro, var):
        """
            Set up some info.
            @param seq: The sequence.
            @param mat: The match data.
            @param nam: The names of the matches.
            @param pro: The protein data for each match.
            @param var: The variant data for each match.
        """
        self.sequence = seq
        """The sequence of the peptide. String."""
        self.matches = mat
        """The raw match data for the peptide. List of Tuples of four ints."""
        self.names = nam
        """The names of each match. List of str."""
        self.protData = pro
        """The protein data for each match. List of list of str."""
        self.varData = var
        """The variant data for each match. List of list of str."""
        # extract stuff from the protein data (start, end, chromo, strand)
        protLocs = []
        for cpd in pro:
            if cpd is None:
                protLocs.append([])
                continue
            ldatS = cpd.index("LOCS{") + 1
            ldatE = cpd.index("}SCOL")
            if (ldatE - ldatS) % 4 != 0:
                raise IOError("Location data must have four entries (chromosome, strand, start, end) per exon.")
            curLocs = []
            for i in range(ldatS, ldatE, 4):
                curLocs.append( (cpd[i], int(cpd[i+1]), int(cpd[i+2]), cpd[i+3] != "-") )
            protLocs.append(curLocs)
        self.protLocs = protLocs
        """The locations of the reference exons of each protein match."""
        # extract stuff from the variant data (cigar, individuals, relevant snps)
        varCigs = ["???" if (cpd is None) else cpd[2] for cpd in var]
        varInds = []
        varDiff = []
        for cpd in var:
            if cpd is None:
                varInds.append([])
                varDiff.append([])
                continue
            idatS = cpd.index("INDS{") + 1
            idatE = cpd.index("}SDNI")
            varInds.append(cpd[idatS:idatE])
            if "DIFS{" in cpd:
                idatS = cpd.index("DIFS{") + 1
                idatE = cpd.index("}SFID")
                relDifs = []
                for cdStr in cpd[idatS:idatE]:
                    relDifs.extend(parseProteosDifferences(cdStr))
                varDiff.append(relDifs)
        self.varCigs = varCigs
        """The cigars for each variant."""
        self.varInds = varInds
        """The individuals having each match."""
        self.varDiff = varDiff
        """The set of differences that could produce each match."""
    def limitToProteinSet(self, protNames):
        """
            Return a copy of this info, limitted to the given proteins.
            @param protNames: The (base) names of the proteins to limit to (set of str).
            @return: The limitted variant of this peptide data.
        """
        subSeq = self.sequence
        subMat = []
        subNam = []
        subPro = []
        subVar = []
        for i in range(len(self.matches)):
            cnam = self.names[i].split("_")[0]
            if not (cnam in protNames):
                continue
            subMat.append(self.matches[i])
            subNam.append(self.names[i])
            subPro.append(self.protData[i])
            subVar.append(self.varData[i])
        return ProteosPeptideInfo(subSeq, subMat, subNam, subPro, subVar)
    def getFullPeptideOwnership(self):
        """
            Get the names of everybody who has this peptide.
            @return: Set of those individuals.
        """
        toRet = set()
        for vinds in self.varInds:
            toRet.update(vinds)
        return toRet
    def getFullPeptideVariants(self):
        """
            Get the set of differences that could produce the peptide.
            @return: Set of those differences.
        """
        toRet = set()
        for vinds in self.varDiff:
            toRet.update(vinds)
        return toRet
    def limitVariantsToNearMatch(self, relDist):
        """
            Limit variants to be within some distance of the peptide match.
            @param relDist: The distance to tolerate.
            @return: A copy or varDiff with the variant data filtered.
        """
        winDiff = []
        for i in range(len(self.matches)):
            pepSVar = self.matches[i][2]
            pepEVar = self.matches[i][3]
            # figure out what is acceptable
            varCig = self.varCigs[i]
            if varCig == "???":
                varCig = repr(pepEVar) + "M"
            cigReg = cigarStringToReferencePositions(0, varCig)
            if (cigReg[0] > 0) or (cigReg[1] > 0):
                raise ValueError("Soft clip in protein cigar.")
            if len(cigReg[2]) < pepEVar:
                raise ValueError("Search reports match beyond end of the cigar.")
            cigReg = cigReg[2]
            lowAcc = None
            higAcc = None
            for j in range(pepSVar, pepEVar):
                if cigReg[j] >= 0:
                    lowAcc = cigReg[j]
                    break
            for j in range(pepEVar-1, pepSVar-1, -1):
                if cigReg[j] >= 0:
                    higAcc = cigReg[j]
                    break
            if lowAcc is None:
                for j in range(pepSVar - 1, -1, -1):
                    if cigReg[j] >= 0:
                        lowAcc = cigReg[j]
                        break
            if higAcc is None:
                for j in range(pepEVar, len(cigReg)):
                    if cigReg[j] >= 0:
                        higAcc = cigReg[j]
                        break
            # get the variants that are acceptable
            winVar = []
            if (not (higAcc is None)) or (not (lowAcc is None)):
                if higAcc is None:
                    higAcc = lowAcc
                if lowAcc is None:
                    lowAcc = higAcc
                lowAcc = lowAcc - relDist
                higAcc = higAcc + relDist
                for cvdat in self.varDiff[i]:
                    if cvdat.protRange is None:
                        continue
                    if cvdat.protRange[1] < lowAcc:
                        continue
                    if higAcc < cvdat.protRange[0]:
                        continue
                    winVar.append(cvdat)
            winDiff.append(winVar)
        return winDiff
    def getAllPeptideGenomeLocations(self):
        """
            Get the genome locations of every peptide match in this entry.
            @return: The set of such entries.
        """
        toRet = []
        for i in range(len(self.matches)):
            toRet.append(self._getPeptideGenomeLocation(i))
        return toRet
    def _getPeptideGenomeLocation(self, index):
        """
            Get the location of a single peptide.
            @param index: The index of said peptide.
            @return: The string representation of said location, and any warnings (empty list if none).
        """
        pepLocGen = ""
        allWarn = []
        matName = self.names[index]
        curMatS = self.matches[index][2]
        curMatE = self.matches[index][3]
        curPLData = self.protLocs[index]
        curPCig = self.varCigs[index]
        # turn the location data into triples (fullLocs)
        if len(curPLData) == 0:
            allWarn.append(matName + " has no location data.")
            return [pepLocGen, allWarn]
        chrom = curPLData[0][0]
        isPos = curPLData[0][3]
        haveWarn = False
        fullLocs = []
        for ldat in curPLData:
            if ldat[3] != isPos:
                allWarn.append(matName + " transcript uses both strands.")
                break
            for i in range(ldat[1], ldat[2]):
                fullLocs.append(i)
                if (ldat[0] != chrom) and not haveWarn:
                    allWarn.append(matName + " uses multiple chromosomes in location data (" + chrom + " and " + ldat[0])
                    haveWarn = True
        if not isPos:
            fullLocs = list(reversed(fullLocs))
        tmpLocs = []
        i = 0
        while i < len(fullLocs):
            if i+3 > len(fullLocs):
                allWarn.append(matName + " location data does not cover a whole number of codons")
                break
            tmpLocs.append([fullLocs[i], fullLocs[i+1], fullLocs[i+2]])
            i = i + 3
        fullLocs = tmpLocs
        # figure out mapping from variant to reference
        tmpCig = [i for i in range(curMatE)]
        if curPCig == "???":
            allWarn.append(matName + " missing cigar")
        else:
            try:
                cigReg = cigarStringToReferencePositions(0, curPCig)
                if (cigReg[0] > 0) or (cigReg[1] > 0):
                    raise ValueError("Soft clip in protein cigar.")
                if len(cigReg[2]) < curMatE:
                    raise ValueError("Search reports match beyond end of the cigar.")
                tmpCig = cigReg[2]
            except ValueError as err:
                errMess = getattr(err, 'message', repr(err))
                allWarn.append(matName + " bad cigar " + curPCig + ": " + errMess)
        # make sure fullLocs has the right length
        haveWarnCigGff = False
        for i in range(len(tmpCig)):
            curCI = tmpCig[i]
            if curCI >= len(fullLocs):
                if not haveWarnCigGff:
                    haveWarnCigGff = True
                    allWarn.append(matName + " cigar extends beyond location data")
                tmpCig[i] = -1
        # peptide_location_genome: chr:start-end  same as before
        pepLocs = []
        for i in range(curMatS,curMatE):
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
        return [pepLocGen, allWarn]


class ProteosPeptideSearch:
    """A search for peptides."""
    def __init__(self, suffArr, lookFor, workFold):
        """
            Set up a search.
            @param suffArr: The suffix array information (ProteosSuffixWrapper).
            @param lookFor: The peptides to search for. Either List of string, or a string (name of a fasta file).
            @param workFold: The folder to store temporary files in.
        """
        lookForFa = None
        if isinstance(lookFor, str):
            lookForFa = lookFor
        else:
            lookForFa = os.path.join(workFold, "searchList.fa")
        searchResFN = os.path.join(workFold, "searchInds.bin")
        nameResFN = os.path.join(workFold, "searchNames.tsv")
        protDataFN = os.path.join(workFold, "protdata.tsv")
        varDataFN = os.path.join(workFold, "vardata.tsv")
        # run search if not already there
        self.needKill = []
        if not os.path.exists(searchResFN):
            self.needKill = [searchResFN, nameResFN, protDataFN, varDataFN]
            if not isinstance(lookFor, str):
                self.needKill.append(lookForFa)
            try:
                runPeptideSearch(suffArr, lookFor, workFold, open(os.devnull, "w"))
            except:
                for subKill in self.needKill:
                    if os.path.exists(subKill):
                        os.remove(subKill)
                raise
        # open all the files
        self.seqFile = open(lookForFa, "r")
        """Get the peptide sequences."""
        self.matFile = open(searchResFN, "rb")
        """Get the matches."""
        self.namFile = open(nameResFN, "r")
        """Get the names of the proteins matched."""
        self.proFile = open(protDataFN, "r")
        """Get the relevant protein data."""
        self.varFile = open(varDataFN, "r")
        """Get the relevant variant file."""
        self.prevSeq = None
        """This will need to read ahead in the fasta: allow saving the lookahead."""
        self.prevMat = None
        """This will need to read ahead in the match file to find things: allow saving the lookahead."""
        self.prevVar = None
        """This may need to read ahead in the variant file to find things: allow saving the lookahead."""
        self.prevPro = None
        """This may need to read ahead in the protein file to find things: allow saving the lookahead."""
        self.curPep = 0
        """The current peptide."""
        self.curMatch = 0
        """The next match index to look for."""
        self.hitSeqEOF = False
        """Whether the end of the peptide file has been hit."""
        self.hitMatEOF = False
        """Whether the end of the match file has been hit."""
        self.hitVarEOF = False
        """Whether the end of the variant file has been hit."""
        self.hitProEOF = False
        """Whether the end of the protein file has been hit."""
    def getNextMatch(self):
        """
            Read the next match.
            @return: The relevant peptide, or null if no such peptide.
        """
        if self.hitSeqEOF:
            return None
        # read the next peptide
        curPSeq = self._readNextPeptideSequence()
        if curPSeq is None:
            return None
        # read matches until no longer right index
        curPMat = self._readRelevantMatches(self.curPep)
        # get the names
        curPNam = self._readSomeLines(self.namFile, len(curPMat), "protein name")
        # get the protein and variant data
        availPro = self._readSomeProteinData(self.curMatch + len(curPMat))
        availVar = self._readSomeVariantData(self.curMatch + len(curPMat))
        curPPro = [None for m in curPMat]
        curPVar = [None for m in curPMat]
        for ai in availPro:
            curPPro[ai - self.curMatch] = availPro[ai]
        for ai in availVar:
            curPVar[ai - self.curMatch] = availVar[ai]
        #package
        toRet = ProteosPeptideInfo(curPSeq, curPMat, curPNam, curPPro, curPVar)
        self.curPep = self.curPep + 1
        self.curMatch = self.curMatch + len(curPMat)
        return toRet
    def _readNextPeptideSequence(self):
        """
            Read the next peptide.
            @return: The relevant sequence.
        """
        while self.prevSeq is None:
            curLin = self.seqFile.readline()
            if curLin == "":
                self.hitSeqEOF = True
                break
            if curLin.strip() == "":
                continue
            self.prevSeq = curLin.strip()
        if self.hitSeqEOF:
            return None
        if self.prevSeq[0] != ">":
            raise IOError("Malformed fasta file.")
        pepSNam = self.prevSeq[1:]
        self.prevSeq = None
        pepSeq = []
        while True:
            curLin = self.seqFile.readline()
            if curLin == "":
                self.hitSeqEOF = True
                break
            curLinS = curLin.strip()
            if curLinS == "":
                continue
            if curLinS[0] == ">":
                self.prevSeq = curLinS
                break
            pepSeq.append(curLinS)
        return "".join(pepSeq)
    def _readRelevantMatches(self, forPInd):
        """
            Read in match data for a given index.
            @param forPInd: The index to read for.
            @return: The relevant match entries.
        """
        allMatch = []
        if self.hitMatEOF:
            return allMatch
        while True:
            curMS = self.prevMat
            if curMS is None:
                curMS = self.matFile.read(32)
                if len(curMS) == 0:
                    self.hitMatEOF = True
                    break
                if len(curMS) != 32:
                    raise IOError("Match file truncated.")
                curMS = struct.unpack(">QQQQ", curMS)
            if curMS[0] < forPInd:
                raise IOError("Match file must be sorted.")
            if curMS[0] > forPInd:
                self.prevMat = curMS
                break
            self.prevMat = None
            allMatch.append(curMS)
        return allMatch
    def _readSomeLines(self, fromStr, numLin, errRep):
        """
            Read some non-empty lines from a file.
            @param fromStr: The file to read from.
            @param numLin: The number of lines.
            @param errRep: The name of the stream, used for errors.
            @return: List of string.
        """
        toRet = []
        while len(toRet) < numLin:
            curLin = fromStr.readline()
            if curLin == "":
                raise IOError("Truncated " + errRep + " file.")
            curLinS = curLin.strip()
            if curLinS == "":
                continue
            toRet.append(curLinS)
        return toRet
    def _readSomeProteinData(self, upToMatchI):
        """
            Read protein match data.
            @param upToMatchI: The match index to read to.
            @return: Map from int to string[].
        """
        allRes = {}
        if self.hitProEOF:
            return allRes
        while True:
            curMS = self.prevPro
            if curMS is None:
                curMS = self.proFile.readline()
                if curMS == "":
                    self.hitProEOF = True
                    break
                if curMS.strip() == "":
                    continue
                curMS = curMS.split()
            if int(curMS[0]) >= upToMatchI:
                self.prevPro = curMS
                break
            self.prevPro = None
            allRes[int(curMS[0])] = curMS
        return allRes
    def _readSomeVariantData(self, upToMatchI):
        """
            Read variant match data.
            @param upToMatchI: The match index to read to.
            @return: Map from int to string[].
        """
        allRes = {}
        if self.hitVarEOF:
            return allRes
        while True:
            curMS = self.prevVar
            if curMS is None:
                curMS = self.varFile.readline()
                if curMS == "":
                    self.hitVarEOF = True
                    break
                if curMS.strip() == "":
                    continue
                curMS = curMS.split()
            if int(curMS[0]) >= upToMatchI:
                self.prevVar = curMS
                break
            self.prevVar = None
            allRes[int(curMS[0])] = curMS
        return allRes
    def closeMatchStream(self):
        """Close down all the data."""
        self.seqFile.close()
        self.matFile.close()
        self.namFile.close()
        self.proFile.close()
        self.varFile.close()
        for killFN in self.needKill:
            os.remove(killFN)
        self.needKill = []

