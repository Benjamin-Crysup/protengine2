#include "profinman_task.h"

#include <string.h>
#include <stdexcept>

#include "whodun_sort.h"
#include "whodun_nmcy.h"
#include "whodun_parse.h"
#include "whodun_datread.h"
#include "whodun_compress.h"
#include "whodun_stringext.h"
#include "whodun_parse_seq.h"

ProfinmanGetMatchRegion::ProfinmanGetMatchRegion(){
	dumpBaseName = 0;
	matchName = 0;
	numPre = 5;
	numPost = 5;
	mySummary = "  Get sequence at/near a match.";
	myMainDoc = "Usage: profinman extfin [OPTION] [FILE]*\n"
		"Get bases neighboring a match.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan extfin 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &dumpBaseName, 0, "    The reference file searched through.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta matchMeta("Match File");
		matchMeta.isFile = true;
		matchMeta.fileExts.insert(".bin");
		addStringOption("--match", &matchName, 0, "    The binary list of match information.\n    --match File.bin\n", &matchMeta);
	ArgumentParserIntMeta preMeta("Prior Bases");
		addIntegerOption("--pre", &numPre, 0, "    How many bases before the match to get.\n    --pre 5\n", &preMeta);
	ArgumentParserIntMeta postMeta("Post Bases");
		addIntegerOption("--post", &numPost, 0, "    How many bases after the match to get.\n    --post 5\n", &postMeta);
	ArgumentParserStrMeta outMeta("Match Region Output File");
		outMeta.isFile = true;
		outMeta.fileWrite = true;
		outMeta.fileExts.insert(".fa");
		outMeta.fileExts.insert(".fasta");
		addStringOption("--out", &outputName, 0, "    The place to write the sequence.\n    --out File.fa\n", &outMeta);
}

int ProfinmanGetMatchRegion::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if((matchName == 0) || (strlen(matchName)==0)){
		matchName = 0;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(numPre < 0){
		argumentError = "Cannot get negative bases.";
		return 1;
	}
	if(numPost < 0){
		argumentError = "Cannot get negative bases.";
		return 1;
	}
	return 0;
}

void ProfinmanGetMatchRegion::runThing(){
	InStream* saveIS = 0;
	OutStream* saveOS = 0;
	SequenceWriter* saveOSS = 0;
	try{
		//get the output file ready
			openSequenceFileWrite(outputName ? outputName : "-", &saveOS, &saveOSS);
		//open up the input
			if(!matchName || (strcmp(matchName,"-")==0)){
				saveIS = new ConsoleInStream();
			}
			else{
				saveIS = new FileInStream(matchName);
			}
		//open up the reference
			std::string baseFN(dumpBaseName);
			std::string blockFN = baseFN + ".blk";
			std::string fastiFN = baseFN + ".fai";
			GZipCompressionMethod compMeth;
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			GailAQSequenceReader gfaIn(&blkComp, fastiFN.c_str());
			uintptr_t numEntries = gfaIn.getNumEntries();
		//place to store a name
			std::string matchEntName;
			std::string matchEntPreName;
			std::string matchEntPostName;
		//start reading matches
			uintptr_t numEnts = 0;
			char entryBuff[MATCH_ENTRY_SIZE];
			uintptr_t numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
			while(numRead){
				//parse the match entry
					if(numRead != MATCH_ENTRY_SIZE){
						throw std::runtime_error("Incomplete match at end of file.");
					}
					//uintptr_t lookFor = be2nat64(entryBuff);
					uintptr_t foundIn = be2nat64(entryBuff+8);
					uintptr_t foundAt = be2nat64(entryBuff+16);
					uintptr_t foundTo = be2nat64(entryBuff+24);
				//idiot checks
					if(foundTo < foundAt){
						throw std::runtime_error("Bad match entry (high index below low index).");
					}
					if(foundIn >= numEntries){
						throw std::runtime_error("Bad match entry (reference sequence index too big).");
					}
					uintptr_t foundInLen = gfaIn.getEntryLength(foundIn);
					if(foundTo > foundInLen){
						throw std::runtime_error("Bad match entry (match extends beyond reference sequence).");
					}
				//get the sequence
					uintptr_t postTo = foundTo + numPost;
						postTo = std::min(postTo, foundInLen);
					uintptr_t preAt = std::max((intptr_t)0, (intptr_t)foundAt - numPre);
					gfaIn.getEntrySubsequence(foundIn, preAt, postTo);
				//make a name
					char numBuff[4*sizeof(uintmax_t)+4];
					sprintf(numBuff, "%ju", (uintmax_t)numEnts);
					matchEntName.clear();
						matchEntName.append("match_");
						matchEntName.append(numBuff);
					matchEntPreName.clear();
						matchEntPreName.append(matchEntName);
						matchEntPreName.append("_pre");
					matchEntPostName.clear();
						matchEntPostName.append(matchEntName);
						matchEntPostName.append("_post");
				//output pre, match and post
					saveOSS->nextNameLen = matchEntPreName.size();
						saveOSS->nextShortNameLen = matchEntPreName.size();
						saveOSS->nextName = &(matchEntPreName[0]);
						saveOSS->nextSeqLen = foundAt - preAt;
						saveOSS->nextSeq = gfaIn.lastReadSeq;
						saveOSS->nextHaveQual = 0;
						saveOSS->writeNextEntry();
					saveOSS->nextNameLen = matchEntName.size();
						saveOSS->nextShortNameLen = matchEntName.size();
						saveOSS->nextName = &(matchEntName[0]);
						saveOSS->nextSeqLen = foundTo - foundAt;
						saveOSS->nextSeq = gfaIn.lastReadSeq + (foundAt - preAt);
						saveOSS->nextHaveQual = 0;
						saveOSS->writeNextEntry();
					saveOSS->nextNameLen = matchEntPostName.size();
						saveOSS->nextShortNameLen = matchEntPostName.size();
						saveOSS->nextName = &(matchEntPostName[0]);
						saveOSS->nextSeqLen = postTo - foundTo;
						saveOSS->nextSeq = gfaIn.lastReadSeq + (foundTo - preAt);
						saveOSS->nextHaveQual = 0;
						saveOSS->writeNextEntry();
				numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
				numEnts++;
			}
	}
	catch(std::exception& err){
		if(saveIS){ delete(saveIS); }
		if(saveOS){ delete(saveOS); }
		if(saveOS){ delete(saveOSS); }
		throw;
	}
	if(saveIS){ delete(saveIS); }
	if(saveOS){ delete(saveOS); }
	if(saveOS){ delete(saveOSS); }
}

/**Compare search results.*/
bool compareBinarySearchData(void* unif, void* itemA, void* itemB){
	return memcmp(itemA, itemB, MATCH_ENTRY_SIZE) < 0;
}

/**Compare the locations of search results.*/
bool compareBinarySearchLocationData(void* unif, void* itemA, void* itemB){
	char* itemAC = (char*)itemA;
	char* itemBC = (char*)itemB;
	return memcmp(itemAC+8, itemBC+8, MATCH_ENTRY_SIZE-8) < 0;
}

ProfinmanSortSearchResults::ProfinmanSortSearchResults(){
	origSRName = 0;
	maxRam = 500000000;
	numThread = 1;
	workFolder = 0;
	outputName = 0;
	mySummary = "  Sort search results..";
	myMainDoc = "Usage: profinman sortfound [OPTION] [FILE]*\n"
		"Sort search results.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan sortfound 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta inMeta("Sort Result File");
		inMeta.isFile = true;
		inMeta.fileExts.insert(".bin");
		addStringOption("--in", &origSRName, 0, "    Specify the data to sort.\n    --in File.bin\n", &inMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    How much ram to use.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserIntMeta threadMeta("Threads");
		addIntegerOption("--thread", &numThread, 0, "    How many threads to use.\n    --thread 1\n", &threadMeta);
	ArgumentParserStrMeta workMeta("Working Folder");
		//TODO add folder to meta
		addStringOption("--work", &workFolder, 0, "    The folder to put temporary files in.\n    --work Folder\n", &workMeta);
	ArgumentParserStrMeta outMeta("Search Result File");
		outMeta.isFile = true;
		outMeta.fileWrite = true;
		outMeta.fileExts.insert(".bin");
		addStringOption("--out", &outputName, 0, "    The place to write the results.\n    --out File.bin\n", &outMeta);
}

ProfinmanSortSearchResults::~ProfinmanSortSearchResults(){}

int ProfinmanSortSearchResults::posteriorCheck(){
	if((origSRName == 0) || (strlen(origSRName)==0)){
		origSRName = 0;
	}
	if((workFolder == 0) || (strlen(workFolder)==0)){
		argumentError = "Need to specify a working folder.";
		return 1;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "Will use at least one byte of ram.";
		return 1;
	}
	if(numThread <= 0){
		argumentError = "Need at least one thread.";
		return 1;
	}
	if(maxRam < 4*MATCH_ENTRY_SIZE){
		maxRam = 4*MATCH_ENTRY_SIZE;
	}
	maxRam = MATCH_ENTRY_SIZE * (maxRam / MATCH_ENTRY_SIZE);
	return 0;
}

void ProfinmanSortSearchResults::runThing(){
	InStream* baseIn = 0;
	OutStream* baseOut = 0;
	try{
		//open the output
		if((outputName == 0) || (strcmp(outputName,"-")==0)){
			baseOut = new ConsoleOutStream();
		}
		else{
			baseOut = new FileOutStream(0, outputName);
		}
		//open the input
		if((origSRName == 0) || (strcmp(origSRName,"-")==0)){
			baseIn = new ConsoleInStream();
		}
		else{
			baseIn = new FileInStream(origSRName);
		}
		//sort
		SortOptions useOpts;
			useOpts.compMeth = compareBinarySearchData;
			useOpts.itemSize = MATCH_ENTRY_SIZE;
			useOpts.maxLoad = maxRam;
			useOpts.numThread = numThread;
			useOpts.useUni = 0;
			useOpts.usePool = 0;
		outOfMemoryMergesort(baseIn, workFolder, baseOut, &useOpts);
	}
	catch(std::exception& err){
		if(baseIn){ delete(baseIn); }
		if(baseOut){ delete(baseOut); }
		throw;
	}
	if(baseIn){ delete(baseIn); }
	if(baseOut){ delete(baseOut); }
}

ProfinmanGetMatchName::ProfinmanGetMatchName(){
	dumpBaseName = 0;
	matchName = 0;
	outputName = 0;
	maxRam = 500000000;
	mySummary = "  Get protein names of matches.";
	myMainDoc = "Usage: profinman nameget [OPTION] [FILE]*\n"
		"Get protein names of matches.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan nameget 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &dumpBaseName, 0, "    The reference file searched through.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta matchMeta("Match File");
		matchMeta.isFile = true;
		addStringOption("--match", &matchName, 0, "    The binary list of match information.\n    --match File.bmatch\n", &matchMeta);
	ArgumentParserStrMeta outMeta("Output File");
		outMeta.isFile = true;
		outMeta.fileWrite = true;
		addStringOption("--out", &outputName, 0, "    The place to write the result.\n    --out File.tsv\n", &outMeta);
}

int ProfinmanGetMatchName::posteriorCheck(){
	if((dumpBaseName==0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if((matchName == 0) || (strlen(matchName)==0)){
		matchName = 0;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "Will use at least one byte of ram.";
		return 1;
	}
	if(maxRam < 4*MATCH_ENTRY_SIZE){
		maxRam = 4*MATCH_ENTRY_SIZE;
	}
	maxRam = MATCH_ENTRY_SIZE * (maxRam / MATCH_ENTRY_SIZE);
	return 0;
}

void ProfinmanGetMatchName::runThing(){
	InStream* saveIS = 0;
	OutStream* saveOS = 0;
	try{
		//open up the input
			if(!matchName || (strcmp(matchName,"-")==0)){
				saveIS = new ConsoleInStream();
			}
			else{
				saveIS = new FileInStream(matchName);
			}
		//open up the output
			if(!outputName || (strcmp(outputName,"-")==0)){
				saveOS = new ConsoleOutStream();
			}
			else{
				saveOS = new FileOutStream(0, outputName);
			}
		//open up the reference
			std::string baseFN(dumpBaseName);
			std::string blockFN = baseFN + ".blk";
			std::string fastiFN = baseFN + ".fai";
			GZipCompressionMethod compMeth;
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			GailAQSequenceReader gfaIn(&blkComp, fastiFN.c_str());
			uintptr_t numEntries = gfaIn.getNumEntries();
		//place to store stuff
			std::vector<char> preloadEnts;
			std::map<uintptr_t, std::pair<uintptr_t,uintptr_t> > entWindices;
			std::string saveFounds;
		//start reading matches
			char entryBuff[MATCH_ENTRY_SIZE];
			uintptr_t numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
			while(numRead || preloadEnts.size()){
				//add the entry to the entities
				if(numRead){
					if(numRead != MATCH_ENTRY_SIZE){
						throw std::runtime_error("Incomplete match at end of file.");
					}
					preloadEnts.insert(preloadEnts.end(), entryBuff, entryBuff + numRead);
				}
				//if enough entries (or nothing left), handle
				if((numRead == 0) || (preloadEnts.size() > (uintptr_t)maxRam)){
					for(uintptr_t i = 0; i<preloadEnts.size(); i+=MATCH_ENTRY_SIZE){
						char* curEnt = &(preloadEnts[i]);
						//uintptr_t lookFor = be2nat64(entryBuff);
						uintptr_t foundIn = be2nat64(curEnt+8);
						//uintptr_t foundAt = be2nat64(curEnt+16);
						//uintptr_t foundTo = be2nat64(curEnt+24);
						//idiot check
							if(foundIn >= numEntries){
								throw std::runtime_error("Bad match entry (reference sequence index too big).");
							}
						//find if not in cache
							if(entWindices.find(foundIn) == entWindices.end()){
								uintptr_t foundInLen = std::min((uintptr_t)1, gfaIn.getEntryLength(foundIn));
								gfaIn.getEntrySubsequence(foundIn, 0, foundInLen);
								entWindices[foundIn] = std::pair<uintptr_t,uintptr_t>(saveFounds.size(), gfaIn.lastReadNameLen);
								saveFounds.insert(saveFounds.end(), gfaIn.lastReadName, gfaIn.lastReadName + gfaIn.lastReadNameLen);
							}
						//get from cache
							std::pair<uintptr_t,uintptr_t> windex = entWindices[foundIn];
							saveOS->writeBytes(&(saveFounds[windex.first]), windex.second);
							saveOS->writeByte('\n');
					}
					preloadEnts.clear();
					entWindices.clear();
					saveFounds.clear();
				}
				//load next
				if(numRead){
					numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
				}
			}
	}
	catch(std::exception& err){
		if(saveIS){ delete(saveIS); }
		if(saveOS){ delete(saveOS); }
		throw;
	}
	if(saveIS){ delete(saveIS); }
	if(saveOS){ delete(saveOS); }
}

#define WHITESPACE " \t\r\n"

/**A rule for a cut.*/
class CuttingRule{
public:
	/**Whether this is a bad rule.*/
	bool badRule;
	/**The length of the start prefix.*/
	int startPreLen;
	/**The valid characters of the start prefix, in reverse order.*/
	char** startPreRev;
	/**Whether all characters must be accounted for in the prefix.*/
	bool startPreAll;
	/**The length of the start suffix.*/
	int startPostLen;
	/**The valid characters of the start suffix.*/
	char** startPost;
	/**The length of the end prefix.*/
	int endPreLen;
	/**The valid characters of the end prefix, in reverse order.*/
	char** endPreRev;
	/**The length of the end suffix.*/
	int endPostLen;
	/**The valid characters of the end suffix.*/
	char** endPost;
	/**Whether all characters must be accounted for in the suffix.*/
	bool endPostAll;
	/**
	 * This will make a cutting rule by parsing a line.
	 * @param parseLine The line to parse.
	 */
	CuttingRule(const char* parseLine){
		//defaults
			badRule = true;
			startPreAll = false;
			endPostAll = false;
			startPreRev = 0;
			startPost = 0;
			endPreRev = 0;
			endPost = 0;
		const char* focChar = parseLine;
		//definitions: simple to use in the presence of errors
			std::vector<char*> allCharClass;
			#define ADVANCE_CHAR focChar++; if(!*focChar){ goto err_handler; }
			#define SKIP_WHITESPACE focChar = focChar + strspn(focChar, WHITESPACE); if(!*focChar){ goto err_handler; }
			#define PARSE_CHAR_CLASS \
				std::string hotChars;\
				char* curAdd;\
				switch(*focChar){\
					case '.':\
						curAdd = (char*)malloc(256);\
						for(int i=1; i<256; i++){curAdd[i-1] = i;}\
						curAdd[255] = 0;\
						allCharClass.push_back(curAdd);\
						break;\
					case '\\':\
						ADVANCE_CHAR\
						curAdd = (char*)malloc(2);\
						curAdd[0] = *focChar;\
						curAdd[1] = 0;\
						allCharClass.push_back(curAdd);\
						break;\
					case '[':\
						ADVANCE_CHAR\
						while(*focChar != ']'){\
							if(*focChar == '\\'){\
								ADVANCE_CHAR\
								hotChars.push_back(*focChar);\
							}\
							else{\
								char hexCA = *focChar;\
								ADVANCE_CHAR\
								char hexCB = *focChar;\
								int hexV = parseHexCode(hexCA, hexCB);\
								if(hexV < 0){ goto err_handler; }\
								hotChars.push_back(hexV);\
							}\
							ADVANCE_CHAR\
						}\
						curAdd = (char*)malloc(hotChars.size() + 1);\
						memcpy(curAdd, hotChars.c_str(), hotChars.size());\
						curAdd[hotChars.size()] = 0;\
						allCharClass.push_back(curAdd);\
						break;\
					case ' ': break;\
					case '\t': break;\
					default:\
						{\
							char hexCA = *focChar;\
							ADVANCE_CHAR\
							char hexCB = *focChar;\
							curAdd = (char*)malloc(2);\
							int hexV = parseHexCode(hexCA, hexCB);\
							if(hexV < 0){ goto err_handler; }\
							curAdd[0] = hexV;\
							curAdd[1] = 0;\
							allCharClass.push_back(curAdd);\
						}\
						break;\
				}\
				ADVANCE_CHAR
			#define BUILD_CHAR_CLASS(lenVar,seqVar,revInd) \
				lenVar = allCharClass.size();\
				seqVar = (char**)malloc(sizeof(char*)*lenVar);\
				for(unsigned ijk = 0; ijk<allCharClass.size(); ijk++){\
					seqVar[revInd] = allCharClass[ijk];\
				}
		SKIP_WHITESPACE
		//see if it needs all
		if(*focChar == '^'){
			startPreAll = true;
			ADVANCE_CHAR
		}
		//get characters until pipe
			while(*focChar != '|'){
				PARSE_CHAR_CLASS
			}
			BUILD_CHAR_CLASS(startPreLen,startPreRev,startPreLen-(ijk+1))
			ADVANCE_CHAR
		//get characters until dash
			allCharClass.clear();
			while(*focChar != '-'){
				PARSE_CHAR_CLASS
			}
			BUILD_CHAR_CLASS(startPostLen,startPost,ijk)
			ADVANCE_CHAR
		//get characters until pipe
			allCharClass.clear();
			while(*focChar != '|'){
				PARSE_CHAR_CLASS
			}
			BUILD_CHAR_CLASS(endPreLen,endPreRev,endPreLen-(ijk+1))
			ADVANCE_CHAR
		//get characters until exclamation or dollar
			allCharClass.clear();
			while(!((*focChar == '!') || (*focChar == '$'))){
				PARSE_CHAR_CLASS
			}
			BUILD_CHAR_CLASS(endPostLen,endPost,ijk)
		//see if the end needs all
		if(*focChar == '$'){
			endPostAll = true;
		}
		badRule = false;
		return;
		
		err_handler:
			for(unsigned ijk=0; ijk<allCharClass.size(); ijk++){free(allCharClass[ijk]);}
	}
	/**Kill the memory.*/
	~CuttingRule(){
		if(startPreRev){
			for(int i = 0; i<startPreLen; i++){
				free(startPreRev[i]);
			}
			free(startPreRev);
		}
		if(startPost){
			for(int i = 0; i<startPostLen; i++){
				free(startPost[i]);
			}
			free(startPost);
		}
		if(endPreRev){
			for(int i = 0; i<endPreLen; i++){
				free(endPreRev[i]);
			}
			free(endPreRev);
		}
		if(endPost){
			for(int i = 0; i<endPostLen; i++){
				free(endPost[i]);
			}
			free(endPost);
		}
	}
};

class CuttingRules{
public:
	/**The rules in question.*/
	std::vector<CuttingRule*> theRules;
	/**
	 * Start a search for matching rules.
	 * @param addLive The place to add live indices to.
	 */
	void startMatchSearch(std::vector<int>* addLive){
		addLive->clear();
		for(unsigned i = 0; i < theRules.size(); i++){
			addLive->push_back(i);
		}
	}
#define SCAN_CODE(lenVar, seqVar, needAGet, forRevIndGet) \
		int splen = strlen(toTest);\
		std::vector<int> nextLive;\
		for(unsigned i = 0; i<toCull->size(); i++){\
			int ci = (*toCull)[i];\
			CuttingRule* curRule = (theRules[ci]);\
			int testLen = curRule->lenVar;\
			char** testVals = curRule->seqVar;\
			bool needAll = needAGet;\
			if(splen < testLen){\
				continue;\
			}\
			if(needAll && (splen > testLen)){\
				continue;\
			}\
			bool foundIt = true;\
			for(int j = 0; j<testLen; j++){\
				char* valChars = testVals[j];\
				int curChar = toTest[forRevIndGet];\
				if(strchr(valChars, curChar) == 0){\
					foundIt = false;\
					break;\
				}\
			}\
			if(foundIt){\
				nextLive.push_back(ci);\
			}\
		}\
		toCull->clear();\
		toCull->insert(toCull->begin(), nextLive.begin(), nextLive.end());
	/**
	 * This will see if the prefix matches the prefix of the starting cut.
	 * @param toCull The live rules to work with.
	 * @param toTest The stuff before the match.
	 */
	void checkStartCutPrefix(std::vector<int>* toCull, const char* toTest){
		SCAN_CODE(startPreLen,startPreRev,curRule->startPreAll,splen-(j+1))
	}
	/**
	 * This will see if the suffix matches the suffix of the starting cut.
	 * @param toCull The live rules to work with.
	 * @param toTest The match.
	 */
	void checkStartCutSuffix(std::vector<int>* toCull, const char* toTest){
		SCAN_CODE(startPostLen,startPost,false,j)
	}
	/**
	 * This will see if the prefix matches the prefix of the ending cut.
	 * @param toCull The live rules to work with.
	 * @param toTest The match.
	 */
	void checkEndCutPrefix(std::vector<int>* toCull, const char* toTest){
		SCAN_CODE(endPreLen,endPreRev,false,splen-(j+1))
	}
	/**
	 * This will see if the suffix matches the suffix of the ending cut.
	 * @param toCull The live rules to work with.
	 * @param toTest The stuff after the match.
	 */
	void checkEndCutSuffix(std::vector<int>* toCull, const char* toTest){
		SCAN_CODE(endPostLen,endPost,curRule->endPostAll,j)
	}
	/**
	 * This will find all cuts the given match satisfies.
	 * @param preMatch The stuff before the match.
	 * @param match The match.
	 * @param postMatch The stuff after the match.
	 * @param toFill The place to put the indices of the matching rules.
	 */
	void getMatchingCuts(const char* preMatch, const char* match, const char* postMatch, std::vector<int>* toFill){
		toFill->clear();
		startMatchSearch(toFill);
		checkStartCutPrefix(toFill, preMatch);
		checkStartCutSuffix(toFill, match);
		checkEndCutPrefix(toFill, match);
		checkEndCutSuffix(toFill, postMatch);
	}
	/**
	 * Loads cutting rules from a file.
	 * @param toParse The file to parse. Closed on return.
	 * @param toRepErr The place to put errors, if any.
	 */
	CuttingRules(std::istream* toParse, std::ostream* toRepErr){
		int lineCount = 0;
		std::string curLine;
		while(std::getline(*toParse, curLine)){
			lineCount++;
			if(strspn(curLine.c_str(), WHITESPACE) == curLine.size()){
				continue;
			}
			CuttingRule* curRule = new CuttingRule(curLine.c_str());
			if(curRule->badRule){
				if(toRepErr){
					*toRepErr << "ERROR: Bad cutting rule on line " << lineCount << ": " << curLine << std::endl;
				}
				delete(curRule);
			}
			else{
				theRules.push_back(curRule);
			}
		}
	}
	/**Kill the memory.*/
	~CuttingRules(){
		for(unsigned i = 0; i<theRules.size(); i++){
			delete (theRules[i]);
		}
	}
};

ProfinmanFilterDigestMatches::ProfinmanFilterDigestMatches(){
	dumpBaseName = 0;
	matchName = 0;
	outputName = 0;
	digestName = 0;
	maxRam = 500000000;
	mySummary = "  Filter matches consistent with a digest.";
	myMainDoc = "Usage: profinman matfildig [OPTION] [FILE]*\n"
		"Filter matches consistent with a digest.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan matfildig 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Reference File");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".gail");
		addStringOption("--ref", &dumpBaseName, 0, "    The reference file searched through.\n    --ref File.gail\n", &dumpMeta);
	ArgumentParserStrMeta matchMeta("Match File");
		matchMeta.isFile = true;
		addStringOption("--match", &matchName, 0, "    The binary list of match information.\n    --match File.bmatch\n", &matchMeta);
	ArgumentParserStrMeta outMeta("Filtered Result File");
		outMeta.isFile = true;
		outMeta.fileWrite = true;
		outMeta.fileExts.insert(".bin");
		addStringOption("--out", &outputName, 0, "    The place to write the results.\n    --out File.bin\n", &outMeta);
	ArgumentParserStrMeta digMeta("Digest Specification");
		digMeta.isFile = true;
		addStringOption("--dig", &digestName, 0, "    The digest specification file.\n    --dig File.dig\n", &digMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    How much ram to use.\n    --ram 500000000\n", &ramMeta);
}

int ProfinmanFilterDigestMatches::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify a reference.";
		return 1;
	}
	if((digestName == 0) || (strlen(digestName)==0)){
		argumentError = "Need to specify a digest to filter on.";
		return 1;
	}
	if((matchName == 0) || (strlen(matchName)==0)){
		matchName = 0;
	}
	if((outputName == 0) || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "Will use at least one byte of ram.";
		return 1;
	}
	if(maxRam < 4*MATCH_ENTRY_SIZE){
		maxRam = 4*MATCH_ENTRY_SIZE;
	}
	maxRam = MATCH_ENTRY_SIZE * (maxRam / MATCH_ENTRY_SIZE);
	return 0;
}

void ProfinmanFilterDigestMatches::runThing(){
	//load in the digest
		std::ifstream ruleStr(digestName);
		if(!ruleStr){
			throw std::runtime_error("Could not open digest file.");
		}
		CuttingRules digestRule(&ruleStr, &(std::cerr));
		ruleStr.close();
	//note the longest prefix and suffix
		uintptr_t numPre = 0;
		uintptr_t numPost = 0;
		for(uintptr_t i = 0; i<digestRule.theRules.size(); i++){
			CuttingRule* curRule = digestRule.theRules[i];
			numPre = std::max(numPre, (uintptr_t)(curRule->startPreLen));
			numPost = std::max(numPost, (uintptr_t)(curRule->endPostLen));
		}
	//start matching
	InStream* saveIS = 0;
	OutStream* saveOS = 0;
	try{
		//open up the output
			saveOS = outputName ? (OutStream*)(new FileOutStream(0, outputName)) : (OutStream*)(new ConsoleOutStream());
		//open up the input
			saveIS = matchName ? (InStream*)(new FileInStream(matchName)) : (InStream*)(new ConsoleInStream());
		//open up the reference
			std::string baseFN(dumpBaseName);
			std::string blockFN = baseFN + ".blk";
			std::string fastiFN = baseFN + ".fai";
			GZipCompressionMethod compMeth;
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			GailAQSequenceReader gfaIn(&blkComp, fastiFN.c_str());
			uintptr_t numEntries = gfaIn.getNumEntries();
		//start reading matches
			std::vector<char> preloadEnts;
			std::string preStore;
			std::string midStore;
			std::string postStore;
			std::vector<int> goodMatch;
			char entryBuff[MATCH_ENTRY_SIZE];
			uintptr_t numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
			while(numRead || preloadEnts.size()){
				//add the entry to the entities
				if(numRead){
					if(numRead != MATCH_ENTRY_SIZE){
						throw std::runtime_error("Incomplete match at end of file.");
					}
					preloadEnts.insert(preloadEnts.end(), entryBuff, entryBuff + numRead);
				}
				//if enough entries (or nothing left), handle
				if((numRead == 0) || (preloadEnts.size() > (uintptr_t)maxRam)){
					//sort them by foundIn, foundAt and foundTo
						SortOptions useOpts;
							useOpts.compMeth = compareBinarySearchLocationData;
							useOpts.itemSize = MATCH_ENTRY_SIZE;
							useOpts.maxLoad = maxRam;
							useOpts.numThread = 1;
							useOpts.useUni = 0;
						inMemoryMergesort(preloadEnts.size() / MATCH_ENTRY_SIZE, &(preloadEnts[0]), &useOpts);
					//run down
					uintptr_t foundInLen = 0;
					uintptr_t lastLoad = -1;
					for(uintptr_t i = 0; i<preloadEnts.size(); i+=MATCH_ENTRY_SIZE){
						char* curEnt = &(preloadEnts[i]);
						//uintptr_t lookFor = be2nat64(entryBuff);
						uintptr_t foundIn = be2nat64(curEnt+8);
						uintptr_t foundAt = be2nat64(curEnt+16);
						uintptr_t foundTo = be2nat64(curEnt+24);
						//idiot checks
							if(foundTo < foundAt){
								throw std::runtime_error("Bad match entry (high index below low index).");
							}
							if(foundIn >= numEntries){
								throw std::runtime_error("Bad match entry (reference sequence index too big).");
							}
							if(lastLoad != foundIn){
								//load if not already in
								lastLoad = foundIn;
								foundInLen = gfaIn.getEntryLength(foundIn);
								gfaIn.getEntrySubsequence(foundIn, 0, foundInLen);
							}
							if(foundTo > foundInLen){
								throw std::runtime_error("Bad match entry (match extends beyond reference sequence).");
							}
						//get the sequence
							preStore.clear(); midStore.clear(); postStore.clear();
							uintptr_t postTo = foundTo + numPost;
								postTo = std::min(postTo, foundInLen);
							uintptr_t preAt = std::max((intptr_t)0, (intptr_t)foundAt - (intptr_t)numPre);
							const char* curSeq = gfaIn.lastReadSeq;
							preStore.insert(preStore.end(), curSeq + preAt, curSeq + foundAt);
							midStore.insert(midStore.end(), curSeq + foundAt, curSeq + foundTo);
							postStore.insert(postStore.end(), curSeq + foundTo, curSeq + postTo);
						//see if it matches
							goodMatch.clear();
							digestRule.getMatchingCuts(preStore.c_str(), midStore.c_str(), postStore.c_str(), &goodMatch);
							if(goodMatch.size()){
								saveOS->writeBytes(curEnt, MATCH_ENTRY_SIZE);
							}
					}
					preloadEnts.clear();
				}
				//load next
				if(numRead){
					numRead = saveIS->readBytes(entryBuff, MATCH_ENTRY_SIZE);
				}
			}
	}
	catch(std::exception& err){
		if(saveIS){ delete(saveIS); }
		if(saveOS){ delete(saveOS); }
		throw;
	}
	if(saveIS){ delete(saveIS); }
	if(saveOS){ delete(saveOS); }
}

