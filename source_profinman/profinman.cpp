
#include <deque>
#include <fstream>
#include <string.h>
#include <algorithm>

#include "whodun_args.h"
#include "whodun_nmcy.h"
#include "whodun_suffix.h"
#include "whodun_datread.h"
#include "whodun_compress.h"
#include "whodun_parse_seq.h"
#include "whodun_stringext.h"

#include "profinman_task.h"

ProfinmanAction::~ProfinmanAction(){}

bool memBlockCompare(const std::pair<const char*,uintptr_t>& itemA, const std::pair<const char*,uintptr_t>& itemB){
	uintptr_t minLen = std::min(itemA.second, itemB.second);
	int cmpRes = memcmp(itemA.first, itemB.first, minLen);
	return (cmpRes < 0) || ((cmpRes == 0) && (itemA.second < itemB.second));
}

bool SingleCharMemBlockCompare::operator ()(const std::pair<const char*,uintptr_t>& itemA, const std::pair<const char*,uintptr_t>& itemB){
	return itemA.first[compInd] < itemB.first[compInd];
}

void outputSearchResult(uintptr_t lookFor, uintptr_t foundIn, uintptr_t foundAt, uintptr_t foundTo, FILE* outputTo, bool asText){
	if(asText){
		fprintf(outputTo, "%ju\t%ju\t%ju\t%ju\n", (uintmax_t)lookFor, (uintmax_t)foundIn, (uintmax_t)foundAt, (uintmax_t)foundTo);
	}
	else{
		char buffOut[MATCH_ENTRY_SIZE];
		nat2be64(lookFor, buffOut);
		nat2be64(foundIn, buffOut+8);
		nat2be64(foundAt, buffOut+16);
		nat2be64(foundTo, buffOut+24);
		fwrite(buffOut, 1, MATCH_ENTRY_SIZE, outputTo);
	}
}

//TODO
//exttest test match region
//unziptab unpack block comp tsv
//sorttab sort block comp tsv

/**
 * Pick out what to run, and run.
 */
int main(int argc, char** argv){
	//make the possible actions
		std::map<std::string,ProfinmanAction*> allActs;
		ProfinmanBlockSequence acts00; allActs["zipfa"] = &acts00;
		ProfinmanDumpSequence acts01; allActs["unzipfa"] = &acts01;
		ProfinmanSearchSequence acts02; allActs["findfa"] = &acts02;
		ProfinmanGetMatchRegion actm00; allActs["extfin"] = &actm00;
		ProfinmanSortSearchResults actm01; allActs["sortfound"] = &actm01;
		ProfinmanGetMatchName actm02; allActs["nameget"] = &actm02;
		ProfinmanFilterDigestMatches actm03; allActs["matfildig"] = &actm03;
		ProfinmanBuildReference actr00; allActs["safa"] = &actr00;
		ProfinmanDumpReference actr01; allActs["dbg_dumpsa"] = &actr01;
		ProfinmanSearchReference actr02; allActs["findsa"] = &actr02;
		ProfinmanPackTable actd00; allActs["ziptab"] = &actd00;
		ProfinmanSortTableCells actd01; allActs["sorttab"] = &actd01;
		ProfinmanSearchSortedTableCells actd02; allActs["findtabs"] = &actd02;
		ProfinmanSearchTableCells actd03; allActs["findtab"] = &actd03;
		ProfinmanSortTableSearchResult actd04; allActs["sortfindtab"] = &actd04;
	//simple help
		if((argc <= 1) || (strcmp(argv[1],"--help")==0) || (strcmp(argv[1],"-h")==0) || (strcmp(argv[1],"/?")==0)){
			std::cout << "Usage: profinman action OPTIONS" << std::endl;
			std::cout << "Finds and manages proteins in a reference." << std::endl;
			std::cout << "Possible actions are:" << std::endl;
			for(std::map<std::string,ProfinmanAction*>::iterator actIt = allActs.begin(); actIt != allActs.end(); actIt++){
				std::cout << actIt->first << std::endl;
				std::cout << actIt->second->mySummary << std::endl;
			}
			return 0;
		}
		if(strcmp(argv[1], "--version")==0){
			std::cout << "ProFinMan 0.0" << std::endl;
			std::cout << "Copyright (C) 2020 UNT Center for Human Identification" << std::endl;
			std::cout << "License LGPLv3+: GNU LGPL version 3 or later\n    <https://www.gnu.org/licenses/lgpl-3.0.html>\nThis is free software: you are free to change and redistribute it.\n";
			std::cout << "There is NO WARRANTY, to the extent permitted by law." << std::endl;
			return 0;
		}
	//pick the action and parse arguments
		if(allActs.find(argv[1]) == allActs.end()){
			std::cerr << "Unknown action " << argv[1] << std::endl;
			return 1;
		}
		ProfinmanAction* theAct = allActs[argv[1]];
		if(theAct->parseArguments(argc-2, argv+2, &std::cout) < 0){
			std::cerr << theAct->argumentError << std::endl;
			return 1;
		}
		if(theAct->needRun == 0){ return 0; }
	//run the stupid thing
		try{
			theAct->runThing();
		}catch(std::exception& err){
			std::cerr << err.what() << std::endl;
			return 1;
		}
	return 0;
}
