#include "profinman_task.h"

#include <deque>
#include <string.h>
#include <stdexcept>
#include <algorithm>

#include "whodun_datread.h"
#include "whodun_compress.h"
#include "whodun_parse_seq.h"

ProfinmanBlockSequence::ProfinmanBlockSequence(){
	dumpBaseName = 0;
	mySummary = "  Prepare protein sequence files for use.";
	myMainDoc = "Usage: profinman zipfa [OPTION] [FILE]*\n"
		"Takes one or more fasta files and builds an index.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan zipfa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--dump", &dumpBaseName, 0, "    Specify the main location to write to.\n    --dump File.gail\n", &dumpMeta);
}

ProfinmanBlockSequence::~ProfinmanBlockSequence(){}

int ProfinmanBlockSequence::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	srcFAs.push_back(argv[0]);
	return 1;
}

void ProfinmanBlockSequence::printExtraGUIInformation(std::ostream* toPrint){
	(*toPrint) << "STRINGVEC\t<|>\tNAME\tInput Files\tFILE\tREAD\t6\t.fa\t.fa.gz\t.fa.gzip\t.fasta\t.fasta.gz\t.fasta.gzip" << std::endl;
}

int ProfinmanBlockSequence::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify an output.";
		return 1;
	}
	if(srcFAs.size() == 0){
		srcFAs.push_back("-");
	}
	return 0;
}

void ProfinmanBlockSequence::runThing(){
	//open up the output
		std::string baseFN(dumpBaseName);
		std::string blockFN = baseFN + ".blk";
		std::string fastiFN = baseFN + ".fai";
		GZipCompressionMethod compMeth;
		BlockCompOutStream blkComp(0, 0x010000, baseFN.c_str(), blockFN.c_str(), &compMeth);
		GailAQSequenceWriter gfaOut(0, &blkComp, fastiFN.c_str());
	//run down the inputs
	for(uintptr_t i = 0; i<srcFAs.size(); i++){
		InStream* saveIS = 0;
		SequenceReader* saveSS = 0;
		openSequenceFileRead(srcFAs[i], &saveIS, &saveSS);
		try{
			while(saveSS->readNextEntry()){
				gfaOut.nextNameLen = saveSS->lastReadNameLen;
				gfaOut.nextShortNameLen = saveSS->lastReadShortNameLen;
				gfaOut.nextName = saveSS->lastReadName;
				gfaOut.nextSeqLen = saveSS->lastReadSeqLen;
				gfaOut.nextSeq = saveSS->lastReadSeq;
				gfaOut.nextHaveQual = saveSS->lastReadHaveQual;
				gfaOut.nextQual = saveSS->lastReadQual;
				gfaOut.writeNextEntry();
				
			}
			if(saveIS){ delete(saveIS); }
			if(saveSS){ delete(saveSS); }
		}
		catch(std::exception& err){
			if(saveIS){ delete(saveIS); }
			if(saveSS){ delete(saveSS); }
			throw;
		}
	}
}


ProfinmanDumpSequence::ProfinmanDumpSequence(){
	dumpBaseName = 0;
	outputName = 0;
	mySummary = "  Dump protein sequence files to fasta.";
	myMainDoc = "Usage: profinman unzipfa [OPTION]\n"
		"Dumps a reference to fasta.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan unzipfa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Dump Reference");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--dump", &dumpBaseName, 0, "    Specify the reference to dump.\n    --dump File.gail\n", &dumpMeta);
	ArgumentParserStrMeta outMeta("Dump Location");
		outMeta.isFile = true;
		outMeta.fileExts.insert(".fa");
		outMeta.fileExts.insert(".fasta");
		addStringOption("--out", &outputName, 0, "    Specify the file to dump to.\n    --out File.fa\n", &outMeta);
}

ProfinmanDumpSequence::~ProfinmanDumpSequence(){}

int ProfinmanDumpSequence::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify the reference to dump.";
		return 1;
	}
	if(outputName && (strlen(outputName)==0)){
		outputName = 0;
	}
	return 0;
}

void ProfinmanDumpSequence::runThing(){
	//open the output
		OutStream* baseOutS = 0;
		SequenceWriter* baseOut = 0;
	try{
		openSequenceFileWrite(outputName ? outputName : "-", &baseOutS, &baseOut);
		//open the input
			std::string baseFN(dumpBaseName);
			std::string blockFN = baseFN + ".blk";
			std::string fastiFN = baseFN + ".fai";
			GZipCompressionMethod compMeth;
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			GailAQSequenceReader gfaIn(&blkComp, fastiFN.c_str());
		//read
			while(gfaIn.readNextEntry()){
				baseOut->nextNameLen = gfaIn.lastReadNameLen;
				baseOut->nextShortNameLen = gfaIn.lastReadShortNameLen;
				baseOut->nextName = gfaIn.lastReadName;
				baseOut->nextSeqLen = gfaIn.lastReadSeqLen;
				baseOut->nextSeq = gfaIn.lastReadSeq;
				baseOut->nextHaveQual = gfaIn.lastReadHaveQual;
				baseOut->nextQual = gfaIn.lastReadQual;
				baseOut->writeNextEntry();
			}
		//close
			delete(baseOut); baseOut = 0;
			delete(baseOutS); baseOutS = 0;
	}
	catch(std::exception& err){
		if(baseOut){ delete(baseOut); }
		if(baseOutS){ delete(baseOutS); }
		throw;
	}
}

ProfinmanSearchSequence::ProfinmanSearchSequence(){
	dumpBaseName = 0;
	searchName = 0;
	outputName = 0;
	maxRam = 500000000;
	txtOut = false;
	mySummary = "  Search for peptides in a sequence file (slow).";
	myMainDoc = "Usage: profinman findfa [OPTION] [FILE]*\n"
		"Takes a fasta file and looks for the entries in a reference.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan findfa 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &dumpBaseName, 0, "    The reference file to search through.\n    --ref File.gail\n", &dumpMeta);
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
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    Specify a target ram usage, in bytes.\n    --ram 500000000\n", &ramMeta);
}

int ProfinmanSearchSequence::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if((searchName == 0) || (strlen(searchName)==0)){
		searchName = 0;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "Need at least one byte of ram.";
		return 1;
	}
	return 0;
}

void ProfinmanSearchSequence::runThing(){
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
	//open up the input
	InStream* saveIS = 0;
	SequenceReader* saveSS = 0;
	uintptr_t totLoadS = 0;
	try{
		openSequenceFileRead(searchName ? searchName : "-", &saveIS, &saveSS);
		std::string allLoadedSeq;
		std::vector<uintptr_t> loadSeqL;
		std::vector< std::pair<const char*,uintptr_t> > compLoadSeq; //used for sorting
		std::map<const char*,uintptr_t> seqLocInd; //go from sequence location to sequence index
		std::deque< std::pair<uintptr_t,uintptr_t> > liveAction; //Current matches
		SingleCharMemBlockCompare charComp; //compare single character
		//read in some amount of sequence
		int moreData = true;
		while(moreData){
			moreData = saveSS->readNextEntry();
			if(moreData){
				loadSeqL.push_back(saveSS->lastReadSeqLen);
				allLoadedSeq.insert(allLoadedSeq.end(), saveSS->lastReadSeq, saveSS->lastReadSeq + saveSS->lastReadSeqLen);
			}
			if((!moreData && allLoadedSeq.size()) || (allLoadedSeq.size() > (uintptr_t)maxRam)){
				//sort the loaded sequences, prepare for main
				compLoadSeq.clear(); seqLocInd.clear();
				uintptr_t curOff = 0;
				for(uintptr_t i = 0; i<loadSeqL.size(); i++){
					uintptr_t curL = loadSeqL[i];
					const char* curStr = &(allLoadedSeq[curOff]);
					compLoadSeq.push_back( std::pair<const char*,uintptr_t>(curStr, curL) );
					seqLocInd[curStr] = totLoadS + i;
					curOff+=curL;
				}
				std::sort(compLoadSeq.begin(), compLoadSeq.end(), memBlockCompare);
				//open up the reference and search
				std::string baseFN(dumpBaseName);
				std::string blockFN = baseFN + ".blk";
				std::string fastiFN = baseFN + ".fai";
				GZipCompressionMethod compMeth;
				BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
				GailAQSequenceReader gfaIn(&blkComp, fastiFN.c_str());
				uintptr_t gailInd = 0;
				while(gfaIn.readNextEntry()){
					for(uintptr_t si = 0; si<gfaIn.lastReadSeqLen; si++){
						//char curC = gfaIn.lastReadSeq[si];
						//add things that start here
							liveAction.push_front( std::pair<uintptr_t,uintptr_t>(0, compLoadSeq.size()) );
						//try to match stuff already on the docket
							uintptr_t curHist = liveAction.size();
							while(curHist){
								curHist--;
								std::pair<uintptr_t,uintptr_t> curPosIs = liveAction[curHist];
								//add matches agains the character
								charComp.compInd = curHist;
								std::pair<const char*,uintptr_t> lowLook(gfaIn.lastReadSeq + (si-curHist), curHist+1);
								std::pair<const char*,uintptr_t> higLook(gfaIn.lastReadSeq + (si-curHist), gfaIn.lastReadSeqLen - (si-curHist));
								std::vector< std::pair<const char*,uintptr_t> >::iterator matchLow = std::lower_bound(compLoadSeq.begin() + curPosIs.first, compLoadSeq.begin() + curPosIs.second, lowLook, charComp);
								std::vector< std::pair<const char*,uintptr_t> >::iterator matchHig = std::upper_bound(compLoadSeq.begin() + curPosIs.first, compLoadSeq.begin() + curPosIs.second, higLook, charComp);
								//remove any matches
								while((matchLow != matchHig) && (matchLow->second == (curHist+1))){
									uintptr_t matSInd = seqLocInd[matchLow->first];
									outputSearchResult(matSInd, gailInd, si-curHist, si+1, dumpTo, txtOut);
									matchLow++;
								}
								//update
								curPosIs.first = matchLow - compLoadSeq.begin();
								curPosIs.second = matchHig - compLoadSeq.begin();
								liveAction[curHist] = curPosIs;
							}
						//if the last thing is empty, kill it
							while(liveAction.size()){
								std::pair<uintptr_t,uintptr_t> curPosIs = liveAction[liveAction.size()-1];
								if(curPosIs.first == curPosIs.second){
									liveAction.pop_back();
								}
								else{
									break;
								}
							}
					}
					liveAction.clear();
					gailInd++;
				}
				//prepare for the next round
				totLoadS += loadSeqL.size();
				allLoadedSeq.clear();
				loadSeqL.clear();
			}
		}
		if(saveIS){ delete(saveIS); }
		if(saveSS){ delete(saveSS); }
		if(killDump){ fclose(dumpTo); }
	}
	catch(std::exception& err){
		if(saveIS){ delete(saveIS); }
		if(saveSS){ delete(saveSS); }
		if(killDump){ fclose(dumpTo); }
		throw;
	}
}


