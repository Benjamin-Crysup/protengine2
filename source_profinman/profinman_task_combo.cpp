#include "profinman_task.h"

#include <iostream>
#include <string.h>
#include <stdexcept>
#include <algorithm>

#include "whodun_sort.h"
#include "whodun_suffix.h"
#include "whodun_datread.h"
#include "whodun_compress.h"
#include "whodun_parse_seq.h"
#include "whodun_stringext.h"

ProfinmanDumpReference::ProfinmanDumpReference(){
	comboName = 0;
	referenceName = 0;
	numPost = 0;
	mySummary = "  Dump a suffix array: debugging.";
	myMainDoc = "Usage: profinman dbg_dumpsa [OPTION]\n"
		"Dump a suffix array.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan dbg_dumpsa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta comboMeta("Combo File");
		comboMeta.isFile = true;
		comboMeta.fileExts.insert(".gail.sa");
		addStringOption("--sa", &comboName, 0, "    The suffix array.\n    --sa File.gail.sa\n", &comboMeta);
	ArgumentParserIntMeta postMeta("Entry Bases");
		addIntegerOption("--base", &numPost, 0, "    How many bases of the entry to output.\n    --base 0\n", &postMeta);
	ArgumentParserStrMeta refMeta("Reference File");
		refMeta.isFile = true;
		refMeta.fileExts.insert(".gail");
		addStringOption("--ref", &referenceName, 0, "    The reference sequence.\n    --ref File.gail\n", &refMeta);
}

ProfinmanDumpReference::~ProfinmanDumpReference(){
}

int ProfinmanDumpReference::posteriorCheck(){
	if((comboName == 0) || (strlen(comboName)==0)){
		argumentError = "Need to specify a suffix array.";
		return 1;
	}
	if(numPost < 0){
		argumentError = "Cannot get a negative number of bases.";
		return 1;
	}
	if((numPost != 0) && ((referenceName == 0) || (strlen(referenceName)==0))){
		argumentError = "If outputing bases, need to specify the reference.";
		return 1;
	}
	return 0;
}

void ProfinmanDumpReference::runThing(){
	const char* hexConv = "0123456789ABCDEF";
	//open up the combo
		std::string cbaseFN(comboName);
		std::string cblockFN = cbaseFN + ".blk";
		GZipCompressionMethod ccompMeth;
		BlockCompInStream comboB(cbaseFN.c_str(), cblockFN.c_str(), &ccompMeth);
		uintptr_t comboEnts = comboB.getUncompressedSize();
		if(comboEnts % COMBO_ENTRY_SIZE){ throw std::runtime_error("Malformed combo file."); }
		comboEnts = comboEnts / COMBO_ENTRY_SIZE;
	//open up the reference, if any
		GZipCompressionMethod rcompMeth;
		BlockCompInStream* rblkComp;
		GailAQSequenceReader* gfaIn;
		if(referenceName){
			std::string rbaseFN(referenceName);
			std::string rblockFN = rbaseFN + ".blk";
			std::string rfastiFN = rbaseFN + ".fai";
			rblkComp = new BlockCompInStream(rbaseFN.c_str(), rblockFN.c_str(), &rcompMeth);
			gfaIn = new GailAQSequenceReader(rblkComp, rfastiFN.c_str());
		}
	//run down the entries
		std::vector<char> postSeq;
		char curEnt[COMBO_ENTRY_SIZE];
		for(uintptr_t i = 0; i<comboEnts; i++){
			comboB.readBytes(curEnt, COMBO_ENTRY_SIZE);
			uintptr_t seqInd = be2nat64(curEnt);
			uintptr_t charIndS = be2nat64(curEnt+8);
			//get the sequence to dump
			postSeq.clear();
			if(numPost){
				uintptr_t seqLen = gfaIn->getEntryLength(seqInd);
				uintptr_t charIndE = charIndS + numPost;
					charIndE = std::min(charIndE, seqLen);
				gfaIn->getEntrySubsequence(seqInd, charIndS, charIndE);
				uintptr_t numGot = charIndE - charIndS;
				for(uintptr_t j = 0; j<numGot; j++){
					char curGot = gfaIn->lastReadSeq[j];
					postSeq.push_back('\t');
					postSeq.push_back(hexConv[(curGot >> 4)&0x0F]);
					postSeq.push_back(hexConv[curGot&0x0F]);
				}
				for(intptr_t j = numGot; j<numPost; j++){
					postSeq.push_back('\t');
					postSeq.push_back('*');
				}
			}
			postSeq.push_back(0);
			std::cout << seqInd << "\t" << charIndS << &(postSeq[0]) << std::endl;
		}
	//clean up
		if(referenceName){
			delete(gfaIn);
			delete(rblkComp);
		}
}

ProfinmanSearchReference::ProfinmanSearchReference(){
	referenceName = 0;
	comboName = 0;
	searchName = 0;
	outputName = 0;
	txtOut = false;
	mySummary = "  Search for peptides in a suffix array.";
	myMainDoc = "Usage: profinman findsa [OPTION] [FILE]*\n"
		"Takes a fasta file and looks for the entries in a reference.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan findsa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &referenceName, 0, "    The reference file to search through.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta comboMeta("Suffix Array File");
		comboMeta.isFile = true;
		comboMeta.fileExts.insert(".gail.sa");
		addStringOption("--sa", &comboName, 0, "    The pre-built suffix array.\n    --sa File.gail.sa\n", &comboMeta);
	ArgumentParserStrMeta entsMeta("Search File");
		entsMeta.isFile = true;
		entsMeta.fileExts.insert(".fasta");
		entsMeta.fileExts.insert(".fa");
		entsMeta.fileExts.insert(".fasta.gzip");
		entsMeta.fileExts.insert(".fa.gzip");
		entsMeta.fileExts.insert(".fasta.gz");
		entsMeta.fileExts.insert(".fa.gz");
		addStringOption("--search", &searchName, 0, "    The file containing sequences to find.\n    --search File.fa\n", &entsMeta);
	ArgumentParserStrMeta outMeta("Search Result File");
		outMeta.isFile = true;
		outMeta.fileWrite = true;
		outMeta.fileExts.insert(".bin");
		addStringOption("--out", &outputName, 0, "    The place to write the results.\n    --out File.bin\n", &outMeta);
	ArgumentParserBoolMeta binMeta("Text Output");
		addBooleanFlag("--text", &txtOut, 1, "    Write out results tsv rather than binary.\n", &binMeta);
}

int ProfinmanSearchReference::posteriorCheck(){
	if(!referenceName || (strlen(referenceName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if(!comboName || (strlen(comboName)==0)){
		argumentError = "Need to specify a suffix array file.";
		return 1;
	}
	if((searchName == 0) || (strlen(searchName)==0)){
		searchName = 0;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	return 0;
}

void ProfinmanSearchReference::runThing(){
	//get the output file ready
	bool killDump = false;
	FILE* dumpTo = 0;
	if((outputName == 0) || (strcmp(outputName, "-")==0)){
		dumpTo = stdout;
	}
	else{
		killDump = true;
		dumpTo = fopen(outputName, "wb");
		if(dumpTo == 0){ throw std::runtime_error("Problem opening output."); }
	}
	//open up the reference
		std::string rbaseFN(referenceName);
		std::string rblockFN = rbaseFN + ".blk";
		std::string rfastiFN = rbaseFN + ".fai";
		GZipCompressionMethod rcompMeth;
		BlockCompInStream rblkComp(rbaseFN.c_str(), rblockFN.c_str(), &rcompMeth);
		GailAQSequenceReader gfaIn(&rblkComp, rfastiFN.c_str());
		uintptr_t numSeqsG = gfaIn.getNumEntries();
	//open up the combo
		std::string cbaseFN(comboName);
		std::string cblockFN = cbaseFN + ".blk";
		GZipCompressionMethod ccompMeth;
		BlockCompInStream comboB(cbaseFN.c_str(), cblockFN.c_str(), &ccompMeth);
		uintptr_t comboEnts = comboB.getUncompressedSize();
		if(comboEnts % COMBO_ENTRY_SIZE){ throw std::runtime_error("Malformed combo file."); }
		comboEnts = comboEnts / COMBO_ENTRY_SIZE;
	//open up the input
		InStream* saveIS = 0;
		SequenceReader* saveSS = 0;
		openSequenceFileRead(searchName ? searchName : "-", &saveIS, &saveSS);
		try{
			uintptr_t curLoadI = 0;
			while(saveSS->readNextEntry()){
				uintptr_t curLen = saveSS->lastReadSeqLen;
				const char* curSeq = saveSS->lastReadSeq;
				//do a lower bound
					uintptr_t lowRangeS = 0;
					uintptr_t countL = comboEnts - lowRangeS;
					while(countL){
						uintptr_t step = countL / 2;
						uintptr_t curTestI = lowRangeS + step;
						//get the entry location
						char entBuff[COMBO_ENTRY_SIZE];
						comboB.seek(COMBO_ENTRY_SIZE*curTestI);
						comboB.readBytes(entBuff, COMBO_ENTRY_SIZE);
						uintptr_t seqInd = be2nat64(entBuff);
						uintptr_t charIndS = be2nat64(entBuff+8);
						uintptr_t charIndE = charIndS + curLen;
						//do the comparison
						if(seqInd >= numSeqsG){ throw std::runtime_error("Suffix array file does not match reference."); }
						bool isLess;
						uintptr_t curSeqLen = gfaIn.getEntryLength(seqInd);
						if(charIndS > curSeqLen){ throw std::runtime_error("Suffix array file does not match reference."); }
						if(charIndE > curSeqLen){
							gfaIn.getEntrySubsequence(seqInd, charIndS, curSeqLen);
							isLess = (memcmp(gfaIn.lastReadSeq, curSeq, (curSeqLen - charIndS)) <= 0);
						}
						else{
							gfaIn.getEntrySubsequence(seqInd, charIndS, charIndE);
							isLess = (memcmp(gfaIn.lastReadSeq, curSeq, curLen) < 0);
						}
						//update the search bounds
						if(isLess){
							lowRangeS = curTestI + 1;
							countL -= (step + 1);
						}
						else{
							countL = step;
						}
					}
				//do an upper bound
					uintptr_t highRangeS = lowRangeS;
					uintptr_t countH = comboEnts - lowRangeS;
					while(countH){
						uintptr_t step = countH / 2;
						uintptr_t curTestI = highRangeS + step;
						//get the entry location
						char entBuff[COMBO_ENTRY_SIZE];
						comboB.seek(COMBO_ENTRY_SIZE*curTestI);
						comboB.readBytes(entBuff, COMBO_ENTRY_SIZE);
						uintptr_t seqInd = be2nat64(entBuff);
						uintptr_t charIndS = be2nat64(entBuff+8);
						uintptr_t charIndE = charIndS + curLen;
						//do the comparison
						if(seqInd >= numSeqsG){ throw std::runtime_error("Suffix array file does not match reference."); }
						bool isLessE;
						uintptr_t curSeqLen = gfaIn.getEntryLength(seqInd);
						if(charIndS > curSeqLen){ throw std::runtime_error("Suffix array file does not match reference."); }
						if(charIndE > curSeqLen){
							gfaIn.getEntrySubsequence(seqInd, charIndS, curSeqLen);
							isLessE = (memcmp(gfaIn.lastReadSeq, curSeq, (curSeqLen - charIndS)) <= 0);
						}
						else{
							gfaIn.getEntrySubsequence(seqInd, charIndS, charIndE);
							isLessE = (memcmp(gfaIn.lastReadSeq, curSeq, curLen) <= 0);
						}
						//update the search bounds
						if(isLessE){
							highRangeS = curTestI + 1;
							countH -= (step + 1);
						}
						else{
							countH = step;
						}
					}
				//report everything in between
					if(highRangeS - lowRangeS){ comboB.seek(COMBO_ENTRY_SIZE*lowRangeS); }
					for(uintptr_t comboI = lowRangeS; comboI < highRangeS; comboI++){
						char entBuff[COMBO_ENTRY_SIZE];
						comboB.readBytes(entBuff, COMBO_ENTRY_SIZE);
						uintptr_t seqInd = be2nat64(entBuff);
						uintptr_t charIndS = be2nat64(entBuff+8);
						outputSearchResult(curLoadI, seqInd, charIndS, charIndS + curLen, dumpTo, txtOut);
					}
				curLoadI++;
			}
		}
		catch(std::exception& err){
			if(saveIS){ delete(saveIS); }
			if(saveSS){ delete(saveSS); }
			if(killDump){ fclose(dumpTo); }
			throw;
		}
		if(saveIS){ delete(saveIS); }
		if(saveSS){ delete(saveSS); }
		if(killDump){ fclose(dumpTo); }
}

ProfinmanMergeReference::ProfinmanMergeReference(){
	refAName = 0;
	refBName = 0;
	comAName = 0;
	comBName = 0;
	referenceName = 0;
	comboName = 0;
	mySummary = "  Merge two suffix arrays (and their gail files).";
	myMainDoc = "Usage: profinman mrgsa [OPTION]\n"
		"Merge two suffix arrays.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan mrgsa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta refAMeta("Reference A");
		refAMeta.isFile = true;
		refAMeta.fileExts.insert(".gail");
		addStringOption("--refA", &refAName, 0, "    The reference for the first suffix array.\n    --refA File.gail\n", &refAMeta);
	ArgumentParserStrMeta refBMeta("Reference B");
		refBMeta.isFile = true;
		refBMeta.fileExts.insert(".gail");
		addStringOption("--refB", &refBName, 0, "    The reference for the second suffix array.\n    --refB File.gail\n", &refBMeta);
	ArgumentParserStrMeta comAMeta("Combo A");
		comAMeta.isFile = true;
		comAMeta.fileExts.insert(".gail.sa");
		addStringOption("--comA", &comAName, 0, "    The first suffix array.\n    --comA File.gail.sa\n", &comAMeta);
	ArgumentParserStrMeta comBMeta("Combo B");
		comBMeta.isFile = true;
		comBMeta.fileExts.insert(".gail");
		addStringOption("--comB", &comBName, 0, "    The second suffix array.\n    --comB File.gail.sa\n", &comBMeta);
	ArgumentParserStrMeta dumpMeta("Reference Out File");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &referenceName, 0, "    The place to write the merged reference.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta comboMeta("Combo Out File");
		comboMeta.isFile = true;
		comboMeta.fileWrite = true;
		comboMeta.fileExts.insert(".gail.sa");
		addStringOption("--out", &comboName, 0, "    The place to write the merged suffix array.\n    --ref File.gail.sa\n", &comboMeta);
}

ProfinmanMergeReference::~ProfinmanMergeReference(){
}

int ProfinmanMergeReference::posteriorCheck(){
	if(!referenceName || (strlen(referenceName)==0)){
		argumentError = "Need to specify an output reference.";
		return 1;
	}
	if(!comboName || (strlen(comboName)==0)){
		argumentError = "Need to specify an output suffix array file.";
		return 1;
	}
	if(!refAName || (strlen(refAName)==0)){
		argumentError = "Need to specify references to merge.";
		return 1;
	}
	if(!refBName || (strlen(refBName)==0)){
		argumentError = "Need to specify references to merge.";
		return 1;
	}
	if(!comAName || (strlen(comAName)==0)){
		argumentError = "Need to specify suffix arrays to merge.";
		return 1;
	}
	if(!comBName || (strlen(comBName)==0)){
		argumentError = "Need to specify suffix arrays to merge.";
		return 1;
	}
	return 0;
}

void ProfinmanMergeReference::runThing(){
	//open the source references
		std::string rbaseAFN(refAName);
		std::string rblockAFN = rbaseAFN + ".blk";
		std::string rfastiAFN = rbaseAFN + ".fai";
		GZipCompressionMethod rcompAMeth;
		BlockCompInStream rblkAComp(rbaseAFN.c_str(), rblockAFN.c_str(), &rcompAMeth);
		GailAQSequenceReader gfaInA(&rblkAComp, rfastiAFN.c_str());
		uintptr_t numSeqsAG = gfaInA.getNumEntries();
		
		std::string rbaseBFN(refBName);
		std::string rblockBFN = rbaseBFN + ".blk";
		std::string rfastiBFN = rbaseBFN + ".fai";
		GZipCompressionMethod rcompBMeth;
		BlockCompInStream rblkBComp(rbaseBFN.c_str(), rblockBFN.c_str(), &rcompBMeth);
		GailAQSequenceReader gfaInB(&rblkBComp, rfastiBFN.c_str());
		//uintptr_t numSeqsBG = gfaInB.getNumEntries();
	//open the source combos
		std::string cbaseAFN(comAName);
		std::string cblockAFN = cbaseAFN + ".blk";
		GZipCompressionMethod ccompAMeth;
		BlockCompInStream comboA(cbaseAFN.c_str(), cblockAFN.c_str(), &ccompAMeth);
		uintptr_t comboAEnts = comboA.getUncompressedSize();
		if(comboAEnts % COMBO_ENTRY_SIZE){ throw std::runtime_error("Malformed combo file."); }
		comboAEnts = comboAEnts / COMBO_ENTRY_SIZE;
		
		std::string cbaseBFN(comBName);
		std::string cblockBFN = cbaseBFN + ".blk";
		GZipCompressionMethod ccompBMeth;
		BlockCompInStream comboB(cbaseBFN.c_str(), cblockBFN.c_str(), &ccompBMeth);
		uintptr_t comboBEnts = comboB.getUncompressedSize();
		if(comboBEnts % COMBO_ENTRY_SIZE){ throw std::runtime_error("Malformed combo file."); }
		comboBEnts = comboBEnts / COMBO_ENTRY_SIZE;
	//open the outputs
		std::string baseFN(referenceName);
		std::string blockFN = baseFN + ".blk";
		std::string fastiFN = baseFN + ".fai";
		GZipCompressionMethod compFAMeth;
		BlockCompOutStream blkCompFA(0, 0x010000, baseFN.c_str(), blockFN.c_str(), &compFAMeth);
		GailAQSequenceWriter gfaOut(0, &blkCompFA, fastiFN.c_str());
		
		GZipCompressionMethod comboCompMeth;
		std::string comFN(comboName);
		std::string comBlkFN = comFN + ".blk";
		BlockCompOutStream blkComp(0, 0x000400, comFN.c_str(), comBlkFN.c_str(), &comboCompMeth);
	//merge the sequences
		while(gfaInA.readNextEntry()){
			gfaOut.nextNameLen = gfaInA.lastReadNameLen;
			gfaOut.nextShortNameLen = gfaInA.lastReadShortNameLen;
			gfaOut.nextName = gfaInA.lastReadName;
			gfaOut.nextSeqLen = gfaInA.lastReadSeqLen;
			gfaOut.nextSeq = gfaInA.lastReadSeq;
			gfaOut.nextHaveQual = gfaInA.lastReadHaveQual;
			gfaOut.nextQual = gfaInA.lastReadQual;
			gfaOut.writeNextEntry();
		}
		while(gfaInB.readNextEntry()){
			gfaOut.nextNameLen = gfaInB.lastReadNameLen;
			gfaOut.nextShortNameLen = gfaInB.lastReadShortNameLen;
			gfaOut.nextName = gfaInB.lastReadName;
			gfaOut.nextSeqLen = gfaInB.lastReadSeqLen;
			gfaOut.nextSeq = gfaInB.lastReadSeq;
			gfaOut.nextHaveQual = gfaInB.lastReadHaveQual;
			gfaOut.nextQual = gfaInB.lastReadQual;
			gfaOut.writeNextEntry();
		}
	//prime the stuff
		char curEntA[COMBO_ENTRY_SIZE];
		char curEntB[COMBO_ENTRY_SIZE];
		uintptr_t numByteA = comboA.readBytes(curEntA, COMBO_ENTRY_SIZE);
		uintptr_t numByteB = comboB.readBytes(curEntB, COMBO_ENTRY_SIZE);
	//load and merge the combos
		uintptr_t lastLoadA = -1;
		uintptr_t lastLoadB = -1;
		while(numByteA && numByteB){
			uintptr_t seqIndA = be2nat64(curEntA);
			uintptr_t charIndSA = be2nat64(curEntA+8);
			uintptr_t seqIndB = be2nat64(curEntB);
			uintptr_t charIndSB = be2nat64(curEntB+8);
			if(seqIndA != lastLoadA){
				uintptr_t readTo = gfaInA.getEntryLength(seqIndA);
				gfaInA.getEntrySubsequence(seqIndA, 0, readTo);
				lastLoadA = seqIndA;
			}
			if(seqIndB != lastLoadB){
				uintptr_t readTo = gfaInB.getEntryLength(seqIndB);
				gfaInB.getEntrySubsequence(seqIndB, 0, readTo);
				lastLoadB = seqIndB;
			}
			if(charIndSA > gfaInA.lastReadSeqLen){ throw std::runtime_error("Suffix array does not match sequence file."); }
			if(charIndSB > gfaInB.lastReadSeqLen){ throw std::runtime_error("Suffix array does not match sequence file."); }
			const char* comSeqA = gfaInA.lastReadSeq + charIndSA;
			const char* comSeqB = gfaInB.lastReadSeq + charIndSB;
			uintptr_t lenSubA = gfaInA.lastReadSeqLen - charIndSA;
			uintptr_t lenSubB = gfaInB.lastReadSeqLen - charIndSB;
			uintptr_t comCompL = std::min(lenSubA, lenSubB);
			int compV = memcmp(comSeqA, comSeqB, comCompL);
			compV = compV ? compV : ((lenSubA < lenSubB) ? -1 : (lenSubA > lenSubB ? 1 : 0));
			if(compV > 0){
				blkComp.writeBytes(curEntA, COMBO_ENTRY_SIZE);
				numByteA = comboA.readBytes(curEntA, COMBO_ENTRY_SIZE);
			}
			else{
				nat2be64(be2nat64(curEntB+8)+numSeqsAG, curEntB+8);
				blkComp.writeBytes(curEntB, COMBO_ENTRY_SIZE);
				numByteB = comboB.readBytes(curEntB, COMBO_ENTRY_SIZE);
			}
		}
		while(numByteA){
			if(numByteA != COMBO_ENTRY_SIZE){
				throw std::runtime_error("Truncated combo file.");
			}
			blkComp.writeBytes(curEntA, COMBO_ENTRY_SIZE);
			numByteA = comboA.readBytes(curEntA, COMBO_ENTRY_SIZE);
		}
		while(numByteB){
			if(numByteB != COMBO_ENTRY_SIZE){
				throw std::runtime_error("Truncated combo file.");
			}
			nat2be64(be2nat64(curEntB+8)+numSeqsAG, curEntB+8);
			blkComp.writeBytes(curEntB, COMBO_ENTRY_SIZE);
			numByteB = comboB.readBytes(curEntB, COMBO_ENTRY_SIZE);
		}
}


//*****************************************************************************
//BUILDING (is a bastard)

/**Thing to pass when sorting.*/
class ProfinmanBuildReferenceSortUni{
public:
	/**The place to input.*/
	InStream* inpFile;
	/**The folder to put temporaries in.*/
	const char* tempFold;
	/**Options for sorting.*/
	SortOptions* sortOpts;
	/**The place to write.*/
	OutStream* outFile;
};
/**Sort stuff.*/
void profinmanBuildReferenceSortTask(void* myUni){
	ProfinmanBuildReferenceSortUni* myUn = (ProfinmanBuildReferenceSortUni*)myUni;
	outOfMemoryMergesort(myUn->inpFile, myUn->tempFold, myUn->outFile, myUn->sortOpts);
}

//***************************
//initial file

/**A task in the initial build.*/
class ProfinmanBuildReferenceInitTask{
public:
	/**The index of this sequence.*/
	uintptr_t seqInd;
	/**The sequence to build for.*/
	std::string forSeq;
};
/**Uniform for the initial make task.*/
class ProfinmanBuildReferenceInitUni{
public:
	/**Get things to turn into output*/
	ThreadProdComCollector<ProfinmanBuildReferenceInitTask>* makeCache;
	/**The place to output.*/
	PreSortMultithreadPipe* dumpFile;
	/**The ID of this task.*/
	uintptr_t taskID;
};
/**Make the initial entries.*/
void profinmanBuildReferenceInitMakeTask(void* myUni){
	std::vector<char> initDump;
	ProfinmanBuildReferenceInitUni* myUn = (ProfinmanBuildReferenceInitUni*)myUni;
	ProfinmanBuildReferenceInitTask* curDo = myUn->makeCache->getThing();
	while(curDo){
		std::string* forSeq = &(curDo->forSeq);
		initDump.clear();
		char curBuildBuff[COMBO_SORT_ENTRY_SIZE];
		//make entries for the main
		nat2be64(curDo->seqInd, curBuildBuff);
		for(uintptr_t i = 1; i<forSeq->size(); i++){
			nat2be64(i-1, curBuildBuff+8);
			nat2be64((*forSeq)[i-1]+1, curBuildBuff+16);
			nat2be64((*forSeq)[i]+1, curBuildBuff+24);
			initDump.insert(initDump.end(), curBuildBuff, curBuildBuff + COMBO_SORT_ENTRY_SIZE);
		}
		//handle the last one special
		nat2be64(forSeq->size()-1, curBuildBuff+8);
		nat2be64((*forSeq)[forSeq->size()-1]+1, curBuildBuff+16);
		nat2be64(0, curBuildBuff+24);
		initDump.insert(initDump.end(), curBuildBuff, curBuildBuff + COMBO_SORT_ENTRY_SIZE);
		//write and return
		myUn->dumpFile->writeBytes(&(initDump[0]), initDump.size());
		myUn->makeCache->taskCache.dealloc(curDo);
		//and continue
		curDo = myUn->makeCache->getThing();
	}
}

//***************************
//rerank

/**Uniform for the rerank task.*/
class ProfinmanBuildReferenceRerankUni{
public:
	/**The offsets. For the first pass, holds the high rank on end. For the second pass, holds the mangle ammount.*/
	uintptr_t useOffset;
	/**The data to edit.*/
	char* toEdit;
	/**The number of entries to edit.*/
	uintptr_t numEdit;
	/**The ID of this task.*/
	uintptr_t taskID;
};
/**First pass rerank.*/
void profinmanBuildReferenceRerankAlpTask(void* myUni){
	ProfinmanBuildReferenceRerankUni* myUn = (ProfinmanBuildReferenceRerankUni*)myUni;
	if(myUn->numEdit == 0){
		myUn->useOffset = 0;
		return;
	}
	//do the first
	char* curFocus = myUn->toEdit;
	uintptr_t prevRank = be2nat64(curFocus + 16);
	uintptr_t prevComp = be2nat64(curFocus + 24);
	uintptr_t nextRank = 0;
	nat2be64(nextRank, curFocus + 16);
	nat2be64(0, curFocus + 24); //compress better
	//do the rest
	uintptr_t numLeft = myUn->numEdit - 1;
	while(numLeft){
		numLeft--;
		curFocus += COMBO_SORT_ENTRY_SIZE;
		uintptr_t curRank = be2nat64(curFocus + 16);
		uintptr_t curComp = be2nat64(curFocus + 24);
		if((curRank != prevRank) || (curComp != prevComp)){
			nextRank++;
			prevRank = curRank;
			prevComp = curComp;
		}
		nat2be64(nextRank, curFocus + 16);
		nat2be64(0, curFocus + 24); //compress better
	}
	myUn->useOffset = nextRank;
}
/**Second pass rerank.*/
void profinmanBuildReferenceRerankBetTask(void* myUni){
	ProfinmanBuildReferenceRerankUni* myUn = (ProfinmanBuildReferenceRerankUni*)myUni;
	if(myUn->numEdit == 0){
		return;
	}
	char* curFocus = myUn->toEdit;
	uintptr_t myOff = myUn->useOffset;
	uintptr_t numDo = myUn->numEdit;
	for(uintptr_t i = 0; i<numDo; i++){
		nat2be64(myOff + be2nat64(curFocus + 16), curFocus+16);
		curFocus += COMBO_SORT_ENTRY_SIZE;
	}
}

//***************************
//next rank

/**A task to get the next rank build.*/
class ProfinmanBuildReferenceNextRankTask{
public:
	/**The bytes to write out.*/
	std::vector<char> initDump;
};
/**Uniform for the initial make task.*/
class ProfinmanBuildReferenceNextRankUni{
public:
	/**The number of zero next ranks at the end.*/
	uintptr_t skipLen;
	/**Get things to turn into output*/
	ThreadProdComCollector<ProfinmanBuildReferenceNextRankTask>* makeCache;
	/**The place to output.*/
	PreSortMultithreadPipe* dumpFile;
	/**The ID of this task.*/
	uintptr_t taskID;
};
/**Fill in the next ranks.*/
void profinmanBuildReferenceNextRankMakeTask(void* myUni){
	ProfinmanBuildReferenceNextRankUni* myUn = (ProfinmanBuildReferenceNextRankUni*)myUni;
	uintptr_t skipLen = myUn->skipLen;
	ProfinmanBuildReferenceNextRankTask* curDo = myUn->makeCache->getThing();
	while(curDo){
		uintptr_t numEnts = curDo->initDump.size() / COMBO_SORT_ENTRY_SIZE;
		if(numEnts > skipLen){
			char* curEnt = &(curDo->initDump[0]);
			numEnts = numEnts - skipLen;
			for(uintptr_t i = 0; i<numEnts; i++){
				nat2be64(be2nat64(curEnt + (COMBO_SORT_ENTRY_SIZE*skipLen) + 16), curEnt + 24);
				curEnt += COMBO_SORT_ENTRY_SIZE;
			}
		}
		//dump the result
		myUn->dumpFile->writeBytes(&(curDo->initDump[0]), curDo->initDump.size());
		myUn->makeCache->taskCache.dealloc(curDo);
		//and continue
		curDo = myUn->makeCache->getThing();
	}
}

//***************************
//showtime

#define THREAD_CACHE_EXTRA 16

ProfinmanBuildReference::ProfinmanBuildReference(){
	referenceName = 0;
	comboName = 0;
	maxRam = 500000000;
	numThread = 1;
	workFolder = 0;
	recoverFile = 0;
	mySummary = "  Build a suffix array of protein sequences.";
	myMainDoc = "Usage: profinman safa [OPTION] [FILE]*\n"
		"Build a suffix array for a sequence file.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan safa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &referenceName, 0, "    The reference to build for.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta comboMeta("Combo Out File");
		comboMeta.isFile = true;
		comboMeta.fileWrite = true;
		comboMeta.fileExts.insert(".gail.sa");
		addStringOption("--out", &comboName, 0, "    The place to write the suffix array.\n    --ref File.gail.sa\n", &comboMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    How much ram to use.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserIntMeta threadMeta("Threads");
		addIntegerOption("--thread", &numThread, 0, "    How many threads to use.\n    --thread 1\n", &threadMeta);
	ArgumentParserStrMeta workMeta("Working Folder");
		workMeta.isFolder = true;
		addStringOption("--work", &workFolder, 0, "    The folder to put temporary files in.\n    --work Folder\n", &workMeta);
	ArgumentParserStrMeta recoMeta("Recover File");
		recoMeta.isFile = true;
		recoMeta.fileWrite = true;
		recoMeta.fileExts.insert(".rec");
		addStringOption("--rec", &recoverFile, 0, "    A recovery file: skip previously finished steps.\n    --rec File.rec\n", &recoMeta);
}

ProfinmanBuildReference::~ProfinmanBuildReference(){
	if(recoverStream){ delete(recoverStream); }
}

int ProfinmanBuildReference::posteriorCheck(){
	if((referenceName == 0) || (strlen(referenceName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if((comboName == 0) || (strlen(comboName)==0)){
		argumentError = "Need to specify the place to write the result.";
		return 1;
	}
	if((workFolder == 0) || (strlen(workFolder)==0)){
		argumentError = "Need to specify a working folder.";
		return 1;
	}
	if(!recoverFile || (strlen(recoverFile)==0)){
		recoverFile = 0;
	}
	if(maxRam <= 0){
		argumentError = "Will use at least one byte of ram.";
		return 1;
	}
	if(numThread <= 0){
		argumentError = "Need at least one thread.";
		return 1;
	}
	if(maxRam < 8*COMBO_SORT_ENTRY_SIZE){
		maxRam = 8*COMBO_SORT_ENTRY_SIZE;
	}
	maxRam = COMBO_SORT_ENTRY_SIZE * (maxRam / COMBO_SORT_ENTRY_SIZE);
	return 0;
}

#define PIPE_BUFFER_SIZE 0x00100000
#define BLOCK_SIZE_INTERNAL 0x010000
#define BLOCK_SIZE_END 0x000400

void ProfinmanBuildReference::runThing(){
	//make threads
		ThreadPool doThreads(numThread);
		ThreadPool prepThread(numThread);
		GZipCompressionMethod baseComp;
	//prepare the sorting methods
		SortOptions rankSortOpts;
			rankSortOpts.itemSize = 4*SUFFIX_ARRAY_CANON_SIZE;
			rankSortOpts.maxLoad = maxRam / 2;
			rankSortOpts.numThread = numThread;
			rankSortOpts.compMeth = MultiStringSuffixRankSortOption_compMeth;
			rankSortOpts.useUni = 0;
			rankSortOpts.usePool = &doThreads;
		SortOptions indSortOpts;
			indSortOpts.itemSize = 4*SUFFIX_ARRAY_CANON_SIZE;
			indSortOpts.maxLoad = maxRam / 2;
			indSortOpts.numThread = numThread;
			indSortOpts.compMeth = MultiStringSuffixIndexSortOption_compMeth;
			indSortOpts.useUni = 0;
			indSortOpts.usePool = &doThreads;
		uintptr_t workRam = maxRam / 2;
		uintptr_t workEntR = COMBO_SORT_ENTRY_SIZE*std::max((uintptr_t)2, workRam / COMBO_SORT_ENTRY_SIZE);
	//load the recovery file
		std::set<std::string> handledTasks;
		if(recoverFile && fileExists(recoverFile)){
			std::ifstream recF(recoverFile);
			std::string curLine;
			while(std::getline(recF, curLine)){
				const char* entStart = curLine.c_str();
					entStart += strspn(entStart, " \t\r\n");
				const char* entEnd = entStart + strcspn(entStart, " \t\r\n");
				if(entStart == entEnd){ continue; }
				handledTasks.insert(std::string(entStart,entEnd));
			}
		}
	//prepare to write the recovery file
		if(recoverFile){ recoverStream = new std::ofstream(recoverFile, std::ofstream::app); }
	//make some names
		std::string workFPrefix = workFolder;
			workFPrefix.append(pathElementSep);
		std::string sortComboChunkA = workFPrefix + "scc_a";
		std::string sortComboChunkB = workFPrefix + "scc_b";
		std::string sortComboChunkAblk = workFPrefix + "scc_a.blk";
		std::string sortComboChunkBblk = workFPrefix + "scc_b.blk";
		std::string sortIndex = workFPrefix + "ssi";
		std::string sortIndexblk = workFPrefix + "ssi.blk";
		std::string refFN(referenceName);
		std::string refBlkFN = refFN + ".blk";
		std::string refFaiFN = refFN + ".fai";
		std::string* curSCC = &sortComboChunkA;
		std::string* nxtSCC = &sortComboChunkB;
		std::string* curSCCB = &sortComboChunkAblk;
		std::string* nxtSCCB = &sortComboChunkBblk;
		std::string* tmpSCC;
		#define SWAP_TARGETS tmpSCC = curSCC; curSCC = nxtSCC; nxtSCC = tmpSCC;    tmpSCC = curSCCB; curSCCB = nxtSCCB; nxtSCCB = tmpSCC;
	//make and sort the initial thing
	uintptr_t totNumString = 0;
	uintptr_t maxStrLen = 0;
	std::vector<uintptr_t> allStrLens;
	{
		if(!(handledTasks.count("init"))){
			{
				MultithreadBlockCompInStream blkComp(refFN.c_str(), refBlkFN.c_str(), &baseComp, numThread, &doThreads);
				GailAQSequenceReader gfaIn(&blkComp, refFaiFN.c_str());
				totNumString = gfaIn.getNumEntries();
				//set up the sort (in its own thread)
					PreSortMultithreadPipe initSPipe(PIPE_BUFFER_SIZE, &doThreads);
					MultithreadBlockCompOutStream initOut(0, BLOCK_SIZE_INTERNAL, nxtSCC->c_str(), nxtSCCB->c_str(), &baseComp, numThread, &doThreads);
					ProfinmanBuildReferenceSortUni initSortU = {&initSPipe, workFolder, &rankSortOpts, &initOut};
					void* sortThread = startThread(profinmanBuildReferenceSortTask, &initSortU);
				//set up the initial producers
					ThreadProdComCollector<ProfinmanBuildReferenceInitTask> makeCache(THREAD_CACHE_EXTRA * numThread);
					std::vector<ProfinmanBuildReferenceInitUni> threadUnis;
					threadUnis.resize(numThread);
					for(intptr_t i = 0; i<numThread; i++){
						ProfinmanBuildReferenceInitUni* curUni = &(threadUnis[i]);
							curUni->makeCache = &makeCache;
							curUni->dumpFile = &initSPipe;
						curUni->taskID = prepThread.addTask(profinmanBuildReferenceInitMakeTask, &(threadUnis[i]));
					}
				//start reading
					#define COMBO_BUILD_INIT_SHUTDOWN \
						makeCache.end();\
						for(intptr_t i = 0; i<numThread; i++){ prepThread.joinTask(threadUnis[i].taskID); }\
						initSPipe.closeWrite();\
						joinThread(sortThread);
					try{
						uintptr_t curStrInd = 0;
						while(gfaIn.readNextEntry()){
							allStrLens.push_back(gfaIn.lastReadSeqLen);
							maxStrLen = std::max(maxStrLen, gfaIn.lastReadSeqLen);
							ProfinmanBuildReferenceInitTask* curTask = makeCache.taskCache.alloc();
							curTask->seqInd = curStrInd;
							curTask->forSeq.clear();
							curTask->forSeq.insert(curTask->forSeq.end(), gfaIn.lastReadSeq, gfaIn.lastReadSeq + gfaIn.lastReadSeqLen);
							makeCache.addThing(curTask);
							curStrInd++;
						}
					}
					catch(std::exception& err){
						COMBO_BUILD_INIT_SHUTDOWN
						throw;
					}
					COMBO_BUILD_INIT_SHUTDOWN
			}
			if(recoverStream){ (*recoverStream) << "init" << std::endl; }
		}
		else{
			GZipCompressionMethod compMeth;
			BlockCompInStream blkComp(refFN.c_str(), refBlkFN.c_str(), &compMeth);
			GailAQSequenceReader gfaIn(&blkComp, refFaiFN.c_str());
			totNumString = gfaIn.getNumEntries();
			//need the total number of strings and the maximum string length
			for(uintptr_t i = 0; i<totNumString; i++){
				uintptr_t curSLen = gfaIn.getEntryLength(i);
				allStrLens.push_back(curSLen);
				maxStrLen = std::max(curSLen, maxStrLen);
			}
		}
	}
	SWAP_TARGETS
	//rerank and resort until everything done
	for(uintptr_t forLen = 4; forLen < 2*maxStrLen; forLen = forLen << 1){
		char numBuff[4*sizeof(uintmax_t)+4];
		sprintf(numBuff, "%ju", forLen);
		std::string curTaskName = "sa_";
			curTaskName.append(numBuff);
		//only do if not previously done
		if(!(handledTasks.count(curTaskName))){
			//basic rerank
			{
				//prepare the sort
				PreSortMultithreadPipe sortindPipe(PIPE_BUFFER_SIZE, &doThreads);
				MultithreadBlockCompOutStream sortindOut(0, BLOCK_SIZE_INTERNAL, sortIndex.c_str(), sortIndexblk.c_str(), &baseComp, numThread, &doThreads);
				ProfinmanBuildReferenceSortUni initSortU = {&sortindPipe, workFolder, &indSortOpts, &sortindOut};
				void* sortThread = startThread(profinmanBuildReferenceSortTask, &initSortU);
				#define COMBO_BUILD_INDSORT_SHUTDOWN \
					sortindPipe.closeWrite();\
					joinThread(sortThread);
				try{
					//rerank things to sort
					uintptr_t chunkPrevRank = -1;
					uintptr_t chunkPrevComp = -1;
					uintptr_t chunkStartOff = 0;
					std::vector<ProfinmanBuildReferenceRerankUni> saveUnis; saveUnis.resize(numThread);
					std::vector<uintptr_t> threadStarts; threadStarts.resize(numThread);
					std::vector<char> curLoad; curLoad.resize(workEntR);
					MultithreadBlockCompInStream initIn(curSCC->c_str(), curSCCB->c_str(), &baseComp, numThread, &doThreads);
					uintptr_t numLoadB = initIn.readBytes(&(curLoad[0]), workEntR);
					while(numLoadB){
						if(numLoadB % COMBO_SORT_ENTRY_SIZE){ throw std::runtime_error("Truncated file."); }
						uintptr_t numLoadE = numLoadB / COMBO_SORT_ENTRY_SIZE;
						//prepare the uniforms
						intptr_t numPerT = numLoadE / numThread;
						intptr_t numExtT = numLoadE % numThread;
						uintptr_t curStart = 0;
						for(intptr_t i = 0; i<numThread; i++){
							ProfinmanBuildReferenceRerankUni* curUni = &(saveUnis[i]);
							curUni->toEdit = &(curLoad[curStart]);
							curUni->numEdit = numPerT + (i<numExtT);
							curStart += (COMBO_SORT_ENTRY_SIZE*curUni->numEdit);
						}
						//note which pieces involve a break
						uintptr_t firstRank = be2nat64(saveUnis[0].toEdit + 16);
						uintptr_t firstComp = be2nat64(saveUnis[0].toEdit + 24);
						threadStarts[0] = (firstRank != chunkPrevRank) || (firstComp != chunkPrevComp);
						chunkPrevRank = be2nat64(&(curLoad[numLoadB-16]));
						chunkPrevComp = be2nat64(&(curLoad[numLoadB-8]));
						for(intptr_t i = 1; i<numThread; i++){
							if(saveUnis[i].numEdit == 0){ threadStarts[i] = 0; continue; }
							ProfinmanBuildReferenceRerankUni* curUni = &(saveUnis[i]);
							if(curUni->numEdit){
								uintptr_t priorRank = be2nat64(curUni->toEdit - 16);
								uintptr_t priorComp = be2nat64(curUni->toEdit - 8);
								uintptr_t curRank = be2nat64(curUni->toEdit + 16);
								uintptr_t curComp = be2nat64(curUni->toEdit + 24);
								threadStarts[i] = ((curRank != priorRank) || (curComp != priorComp));
							}
							else{
								threadStarts[i] = 0;
							}
						}
						//raw-rank in pieces
						for(intptr_t i = 0; i<numThread; i++){
							ProfinmanBuildReferenceRerankUni* curUni = &(saveUnis[i]);
							curUni->taskID = prepThread.addTask(profinmanBuildReferenceRerankAlpTask, curUni);
						}
						for(intptr_t i = 0; i<numThread; i++){ prepThread.joinTask(saveUnis[i].taskID); }
						//adjust each piece
						uintptr_t prevOff = chunkStartOff;
						for(intptr_t i = 0; i<numThread; i++){
							uintptr_t curOffset = prevOff + threadStarts[i];
							prevOff = curOffset + saveUnis[i].useOffset;
							saveUnis[i].useOffset = curOffset;
							saveUnis[i].taskID = prepThread.addTask(profinmanBuildReferenceRerankBetTask, &(saveUnis[i]));
						}
						for(intptr_t i = 0; i<numThread; i++){ prepThread.joinTask(saveUnis[i].taskID); }
						//write
						sortindPipe.writeBytes(&(curLoad[0]), numLoadB);
						//prepare offset rank stuff
						chunkStartOff = prevOff;
						//read next
						numLoadB = initIn.readBytes(&(curLoad[0]), workEntR);
					}
				}
				catch(std::exception& err){
					COMBO_BUILD_INDSORT_SHUTDOWN
					throw;
				}
				COMBO_BUILD_INDSORT_SHUTDOWN
			}
			//make next rank
			{
				//set up the sort
					PreSortMultithreadPipe rrankSPipe(PIPE_BUFFER_SIZE, &doThreads);
					MultithreadBlockCompOutStream rrankOut(0, BLOCK_SIZE_INTERNAL, nxtSCC->c_str(), nxtSCCB->c_str(), &baseComp, numThread, &doThreads);
					ProfinmanBuildReferenceSortUni rrankSortU = {&rrankSPipe, workFolder, &rankSortOpts, &rrankOut};
					void* sortThread = startThread(profinmanBuildReferenceSortTask, &rrankSortU);
				//set up the rerank threads
					MultithreadBlockCompInStream initIn(sortIndex.c_str(), sortIndexblk.c_str(), &baseComp, numThread, &doThreads);
					uintptr_t skipLen = forLen >> 1;
					ThreadProdComCollector<ProfinmanBuildReferenceNextRankTask> makeCache(THREAD_CACHE_EXTRA * numThread);
					std::vector<ProfinmanBuildReferenceNextRankUni> threadUnis;
					threadUnis.resize(numThread);
					for(intptr_t i = 0; i<numThread; i++){
						threadUnis[i] = {skipLen, &makeCache, &rrankSPipe};
						threadUnis[i].taskID = prepThread.addTask(profinmanBuildReferenceNextRankMakeTask, &(threadUnis[i]));
					}
					#define COMBO_BUILD_RERANK_SHUTDOWN \
						makeCache.end();\
						for(intptr_t i = 0; i<numThread; i++){ prepThread.joinTask(threadUnis[i].taskID); }\
						rrankSPipe.closeWrite();\
						joinThread(sortThread);
				//start reading
					try{
						for(uintptr_t i = 0; i<totNumString; i++){
							uintptr_t curStrLen = allStrLens[i];
							ProfinmanBuildReferenceNextRankTask* curTask = makeCache.taskCache.alloc();
							curTask->initDump.resize(COMBO_SORT_ENTRY_SIZE*curStrLen);
							initIn.readBytes(&(curTask->initDump[0]), COMBO_SORT_ENTRY_SIZE*curStrLen);
							makeCache.addThing(curTask);
						}
					}
					catch(std::exception& err){
						COMBO_BUILD_RERANK_SHUTDOWN
						throw;
					}
					COMBO_BUILD_RERANK_SHUTDOWN
			}
			if(recoverStream){ (*recoverStream) << curTaskName << std::endl; }
		}
		SWAP_TARGETS
	}
	//dump to target
	if(!(handledTasks.count("dump"))){
		MultithreadBlockCompInStream endIn(curSCC->c_str(), curSCCB->c_str(), &baseComp, numThread, &doThreads);
		GZipCompressionMethod compMeth;
		std::string comFN(comboName);
		std::string comBlkFN = comFN + ".blk";
		MultithreadBlockCompOutStream blkComp(0, BLOCK_SIZE_END, comFN.c_str(), comBlkFN.c_str(), &baseComp, numThread, &doThreads);
		char curEntBuff[COMBO_SORT_ENTRY_SIZE];
		uintptr_t numR = endIn.readBytes(curEntBuff, COMBO_SORT_ENTRY_SIZE);
		while(numR){
			if(numR != COMBO_SORT_ENTRY_SIZE){ throw std::runtime_error("File truncated."); }
			blkComp.writeBytes(curEntBuff, 16);
			numR = endIn.readBytes(curEntBuff, COMBO_SORT_ENTRY_SIZE);
		}
		if(recoverStream){ (*recoverStream) << "dump" << std::endl; }
	}
	//clean up after yourself
	if(fileExists(sortComboChunkA.c_str())){ killFile(sortComboChunkA.c_str()); }
	if(fileExists(sortComboChunkB.c_str())){ killFile(sortComboChunkB.c_str()); }
	if(fileExists(sortIndex.c_str())){ killFile(sortIndex.c_str()); }
	if(fileExists(sortComboChunkAblk.c_str())){ killFile(sortComboChunkAblk.c_str()); }
	if(fileExists(sortComboChunkBblk.c_str())){ killFile(sortComboChunkBblk.c_str()); }
	if(fileExists(sortIndexblk.c_str())){ killFile(sortIndexblk.c_str()); }
}



