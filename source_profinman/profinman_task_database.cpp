#include "profinman_task.h"

#include <iostream>
#include <string.h>
#include <stdexcept>

#include "whodun_sort.h"
#include "whodun_oshook.h"
#include "whodun_datread.h"
#include "whodun_compress.h"
#include "whodun_parse_table.h"

ProfinmanPackTable::ProfinmanPackTable(){
	dumpBaseName = 0;
	stdoutName[0] = '-'; stdoutName[1] = 0;
	mySummary = "  Pack/compress a tsv database.";
	myMainDoc = "Usage: profinman ziptab [OPTION] [FILE]*\n"
		"Pack/compress a tsv database.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan ziptab 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileExts.insert(".bctsv");
		addStringOption("--dump", &dumpBaseName, 0, "    Specify the main location to write to.\n    --dump File.bctsv\n", &dumpMeta);
}

int ProfinmanPackTable::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	srcFAs.push_back(argv[0]);
	return 1;
}

void ProfinmanPackTable::printExtraGUIInformation(std::ostream* toPrint){
	(*toPrint) << "STRINGVEC\t<|>\tNAME\tInput Files\tFILE\tREAD\t1\t.tsv" << std::endl;
}

ProfinmanPackTable::~ProfinmanPackTable(){}

int ProfinmanPackTable::posteriorCheck(){
	if((dumpBaseName == 0) || (strlen(dumpBaseName)==0)){
		argumentError = "Need to specify an output location.";
		return 1;
	}
	if(srcFAs.size() == 0){
		srcFAs.push_back(stdoutName);
	}
	return 0;
}

void ProfinmanPackTable::runThing(){
	//open the output
		std::string baseFN(dumpBaseName);
		std::string blockFN = baseFN + ".blk";
		std::string fastiFN = baseFN + ".tai";
		GZipCompressionMethod compMeth;
		BlockCompOutStream blkComp(0, 0x010000, baseFN.c_str(), blockFN.c_str(), &compMeth);
		BCompTabularWriter gfaOut(0, &blkComp, fastiFN.c_str());
	//open the TSV
		for(uintptr_t i = 0; i<srcFAs.size(); i++){
			InStream* readF = 0;
			if(strcmp(srcFAs[i],"-")==0){
				readF = new ConsoleInStream();
			}
			else{
				readF = new FileInStream(srcFAs[i]);
			}
			try{
				TSVTabularReader curInT(1, readF);
				while(curInT.readNextEntry()){
					gfaOut.numEntries = curInT.numEntries;
					gfaOut.entrySizes = curInT.entrySizes;
					gfaOut.curEntries = curInT.curEntries;
					gfaOut.writeNextEntry();
				}
				delete(readF); readF = 0;
			}
			catch(std::exception& err){
				if(readF){ delete(readF); }
				throw;
			}
		}
}

ProfinmanSortTableCells::ProfinmanSortTableCells(){
	lookTable = 0;
	outputName = 0;
	maxRam = 500000000;
	numThread = 1;
	stdoutName[0] = '-'; stdoutName[1] = 0;
	mySummary = "  Sort a packed database.";
	myMainDoc = "Usage: profinman sorttab [OPTION] [FILE]*\n"
		"Sort a database.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan sorttab 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".bctsv");
		addStringOption("--dump", &outputName, 0, "    Specify the main location to write to.\n    --dump File.bctsv\n", &dumpMeta);
	ArgumentParserStrMeta origMeta("Unsorted Table");
		origMeta.isFile = true;
		origMeta.fileExts.insert(".bctsv");
		addStringOption("--in", &lookTable, 0, "    Specify the table to sort.\n    --in File.bctsv\n", &origMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    Specify a target ram usage, in bytes.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserIntMeta threadMeta("Threads");
		addIntegerOption("--thread", &numThread, 0, "    How many threads to use.\n    --thread 1\n", &threadMeta);
	ArgumentParserStrVecMeta ocolMeta("Column Codes");
		addStringVectorOption("--col", &indexCols, 0, "    The column indices to sort on, and how to treat them.\n    Prefix with an i for integers, an f for floats\n    --col i0\n", &ocolMeta);
	ArgumentParserStrMeta workMeta("Working Folder");
		workMeta.isFolder = true;
		addStringOption("--work", &workFolder, 0, "    The folder to put temporary files in.\n    --work Folder\n", &workMeta);
}

ProfinmanSortTableCells::~ProfinmanSortTableCells(){
}

/**
 * Check a column specification.
 * @param forSpec The specification.
 * @return An error message.
 */
const char* profdatatsvCheckColumnSpec(std::vector<char*>* forSpec){
	if(forSpec->size() == 0){
		return "Must specify a column.";
	}
	for(uintptr_t i = 0; i<forSpec->size(); i++){
		char* curSpec = (*forSpec)[i];
		if(strlen(curSpec) == 0){
			return "Empty column specification.";
		}
		if(strchr("0123456789",*curSpec)){
			continue;
		}
		if(!strchr("if", *curSpec)){
			return "Unknown type code in column specification.";
		}
		if(strlen(curSpec) == 1){
			return "Empty column specification.";
		}
		if(!strchr("0123456789",curSpec[1])){
			return "Malformed column specification.";
		}
	}
	return 0;
}

/**
 * Check consistency of column specifications.
 * @param forSpecA The specification.
 * @param forSpecB The specification.
 * @return An error message.
 */
const char* profdatatsvCheckConsistent(std::vector<char*>* forSpecA, std::vector<char*>* forSpecB){
	if(forSpecA->size() != forSpecB->size()){
		return "Specifications have different numbers of columns.";
	}
	const char* allowTps = "if";
	for(uintptr_t i = 0; i<forSpecA->size(); i++){
		char* curSpecA = (*forSpecA)[i];
		char* curSpecB = (*forSpecB)[i];
		if(strchr(allowTps, *curSpecA) != strchr(allowTps, *curSpecB)){
			return "Mismatched specification types.";
		}
	}
	return 0;
}

#define COLSPEC_STR 1
#define COLSPEC_INT 2
#define COLSPEC_FLT 3

/**
 * Compile a column specification.
 * @param forSpec THe specification.
 * @param fillInds The indices of the specification.
 * @param fillCode THe way to handle each index.
 * @return Whether there was an error.
 */
const char* profdatatsvCompileColumnSpec(std::vector<char*>* forSpec, std::vector<intptr_t>* fillInds, std::vector<int>* fillCode){
	for(uintptr_t i = 0; i<forSpec->size(); i++){
		char* curFoc = (*forSpec)[i];
		if(*curFoc == 'i'){
			fillCode->push_back(COLSPEC_INT);
			curFoc++;
		}
		else if(*curFoc == 'f'){
			fillCode->push_back(COLSPEC_FLT);
			curFoc++;
		}
		else{
			fillCode->push_back(COLSPEC_STR);
		}
		intptr_t colInd = atol(curFoc);
		if(colInd < 0){
			return "Negative column index.";
		}
		fillInds->push_back(colInd);
	}
	return 0;
}

int ProfinmanSortTableCells::posteriorCheck(){
	if(!lookTable || (strlen(lookTable)==0)){
		argumentError = "Need to specify the table to sort.";
		return 1;
	}
	if(!outputName || (strlen(outputName)==0)){
		argumentError = "Need to specify the output file.";
		return 1;
	}
	if(maxRam <= 0){
		argumentError = "RAM must be positive.";
		return 1;
	}
	if(numThread <= 0){
		argumentError = "Thread count must be positive.";
		return 1;
	}
	if((workFolder == 0) || (strlen(workFolder)==0)){
		argumentError = "Need to specify a working folder.";
		return 1;
	}
	const char* cindErr;
	cindErr = profdatatsvCheckColumnSpec(&indexCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCompileColumnSpec(&indexCols, &fillInds, &fillCode); if(cindErr){ argumentError = cindErr; return 1; }
	return 0;
}

/**Thing to pass when sorting.*/
class ProfinmanSortTablePiecesUniform{
public:
	/**The place to input.*/
	InStream* initSPipe;
	/**The folder to put temporaries in.*/
	const char* workFolder;
	/**Options for sorting.*/
	SortOptions* sortOpts;
	/**The place to write.*/
	OutStream* initOut;
	/**The types of each field.*/
	std::vector<int>* fillCode;
	/**The sizes of each field.*/
	std::vector<uintptr_t>* fieldSizes;
};

/**Run sort from its own thread.*/
void profinmanSortTableThreadFun(void* myUni){
	ProfinmanSortTablePiecesUniform* myU = (ProfinmanSortTablePiecesUniform*)myUni;
	outOfMemoryMergesort(myU->initSPipe, myU->workFolder, myU->initOut, myU->sortOpts);
}

/**Compare mangled table row data.*/
bool profinmanSortTableCompareFun(void* unif, void* itemA,void* itemB){
	char* curFocA = (char*)itemA; curFocA += sizeof(uintptr_t);
	char* curFocB = (char*)itemB; curFocB += sizeof(uintptr_t);
	ProfinmanSortTablePiecesUniform* myU = (ProfinmanSortTablePiecesUniform*)unif;
	for(uintptr_t i = 0; i<myU->fillCode->size(); i++){
		uintptr_t curFSize = (*(myU->fieldSizes))[i];
		int curC = (*(myU->fillCode))[i];
		if(curC == COLSPEC_FLT){
			union {
				double asD;
				char asCs[sizeof(double)];
			} punmeharder;
			memcpy(punmeharder.asCs, curFocA, sizeof(double));
			double cvalA = punmeharder.asD;
			memcpy(punmeharder.asCs, curFocB, sizeof(double));
			double cvalB = punmeharder.asD;
			if(cvalA < cvalB){ return true; }
			if(cvalB < cvalA){ return false; }
		}
		else if(curC == COLSPEC_INT){
			union {
				intptr_t asI;
				char asCs[sizeof(intptr_t)];
			} punmeharder;
			memcpy(punmeharder.asCs, curFocA, sizeof(intptr_t));
			intptr_t cvalA = punmeharder.asI;
			memcpy(punmeharder.asCs, curFocB, sizeof(intptr_t));
			intptr_t cvalB = punmeharder.asI;
			if(cvalA < cvalB){ return true; }
			if(cvalB < cvalA){ return false; }
		}
		else{
			int memcres = memcmp(curFocA, curFocB, curFSize);
			if(memcres < 0){ return true; }
			if(memcres > 0){ return false; }
		}
		curFocA += curFSize;
		curFocB += curFSize;
	}
	return false;
}

#define PIPE_BUFFER_SIZE 4096

/**
 * Packages a table row by index specification.
 * @param entryInd The index of this entry.
 * @param packEnt The row to package.
 * @param packIn The place to store it.
 * @param fromInds The columns of interest.
 * @param fromCodes The types to treat those columns as.
 * @param fieldSizes The number of bytes to use for each column.
 */
void profinmanSortPackageTableRow(uintptr_t entryInd, TabularReader* packEnt, std::vector<char>* packIn, std::vector<intptr_t>* fromInds, std::vector<int>* fromCodes, std::vector<uintptr_t>* fieldSizes){
	{
		union {
			uintptr_t asI;
			char asCs[sizeof(uintptr_t)];
		} punmeharder;
		punmeharder.asI = entryInd;
		packIn->insert(packIn->end(), punmeharder.asCs, punmeharder.asCs + sizeof(uintptr_t));
	}
	for(uintptr_t i = 0; i<fromCodes->size(); i++){
		int curCode = (*fromCodes)[i];
		intptr_t curInd = (*fromInds)[i];
		uintptr_t curFieldS = (*fieldSizes)[i];
		const char* curEnt = packEnt->curEntries[curInd];
		uintptr_t curEntS = packEnt->entrySizes[curInd];
		uintptr_t origSize = packIn->size();
		packIn->insert(packIn->end(), curEnt, curEnt + curEntS);
		if(curCode == COLSPEC_INT){
			packIn->push_back(0);
			union {
				intptr_t asI;
				char asCs[sizeof(intptr_t)];
			} punmeharder;
			punmeharder.asI = atol(&((*packIn)[origSize]));
			packIn->resize(origSize);
			packIn->insert(packIn->end(), punmeharder.asCs, punmeharder.asCs + sizeof(intptr_t));
		}
		else if(curCode == COLSPEC_FLT){
			packIn->push_back(0);
			union {
				double asI;
				char asCs[sizeof(double)];
			} punmeharder;
			punmeharder.asI = atof(&((*packIn)[origSize]));
			packIn->resize(origSize);
			packIn->insert(packIn->end(), punmeharder.asCs, punmeharder.asCs + sizeof(double));
		}
		else{
			if(curEntS > curFieldS){
				packIn->resize(origSize + curFieldS);
			}
			else{
				packIn->insert(packIn->end(), curFieldS - curEntS, 0);
			}
		}
	}
}

/**
 * Repackage entries for a change in size.
 * @param packIn The original thing.
 * @param packOut The repackage.
 * @param fromCodes The type codes for each entry.
 * @param origSizes The original field sizes.
 * @param newSizes The new field sizes.
 */
void profinmanSortRepackageEntries(std::vector<char>* packIn, std::vector<char>* packOut, std::vector<int>* fromCodes, std::vector<uintptr_t>* origSizes, std::vector<uintptr_t>* newSizes){
	if(packIn->size() == 0){ return; }
	char* curFoc = &((*packIn)[0]);
	char* endOrig = curFoc + packIn->size();
	while(curFoc != endOrig){
		//move the index
		packOut->insert(packOut->end(), curFoc, curFoc + sizeof(uintptr_t)); curFoc += sizeof(uintptr_t);
		//move the entries
		for(uintptr_t i = 0; i<fromCodes->size(); i++){
			int curCode = (*fromCodes)[i];
			if(curCode == COLSPEC_INT){
				packOut->insert(packOut->end(), curFoc, curFoc + sizeof(intptr_t)); curFoc += sizeof(intptr_t);
			}
			else if(curCode == COLSPEC_FLT){
				packOut->insert(packOut->end(), curFoc, curFoc + sizeof(double)); curFoc += sizeof(double);
			}
			else{
				uintptr_t origS = (*origSizes)[i];
				uintptr_t newS = (*newSizes)[i];
				if(origS < newS){
					packOut->insert(packOut->end(), curFoc, curFoc + origS); curFoc += origS;
					packOut->insert(packOut->end(), newS - origS, 0);
				}
				else{
					packOut->insert(packOut->end(), curFoc, curFoc + newS); curFoc += origS;
				}
			}
		}
	}
}

void ProfinmanSortTableCells::runThing(){
	std::string baseFN(lookTable);
	std::string blockFN = baseFN + ".blk";
	std::string fastiFN = baseFN + ".tai";
	GZipCompressionMethod compMeth;
	//figure the dumb sizes of each field
		uintptr_t maxCInd = 0;
		std::vector<uintptr_t> fieldSizes;
		fieldSizes.insert(fieldSizes.end(), indexCols.size(), 0);
		for(uintptr_t i = 0; i<indexCols.size(); i++){
			switch(fillCode[i]){
				case COLSPEC_INT:
					fieldSizes[i] = sizeof(intptr_t);
					break;
				case COLSPEC_FLT:
					fieldSizes[i] = sizeof(double);
				default:
					fieldSizes[i] = 0;
			}
			maxCInd = std::max(maxCInd, (uintptr_t)(fillInds[i]));
		}
	//get the largest strings
		{
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			BCompTabularReader gfaOut(&blkComp, fastiFN.c_str());
			while(gfaOut.readNextEntry()){
				if(gfaOut.numEntries <= maxCInd){
					throw std::runtime_error("Table entry has too few columns.");
				}
				for(uintptr_t i = 0; i<indexCols.size(); i++){
					if(fillCode[i] != COLSPEC_STR){ continue; }
					uintptr_t curLen = gfaOut.entrySizes[fillInds[i]];
					fieldSizes[i] = std::max(fieldSizes[i], curLen);
				}
			}
		}
		uintptr_t totItemSize = sizeof(uintptr_t);
		for(uintptr_t i = 0; i<fieldSizes.size(); i++){ totItemSize += fieldSizes[i]; }
	//sort the columns (carrying the indices along for the ride)
		std::string sortIndOutName = workFolder;
			sortIndOutName.append(pathElementSep);
			sortIndOutName.append("stab_ind");
		{
			ProfinmanSortTablePiecesUniform initSortU;
				initSortU.workFolder = workFolder;
				initSortU.fieldSizes = &fieldSizes;
				initSortU.fillCode = &fillCode;
			SortOptions sortOpts;
				sortOpts.itemSize = totItemSize;
				sortOpts.maxLoad = maxRam;
				sortOpts.numThread = numThread;
				sortOpts.compMeth = profinmanSortTableCompareFun;
				sortOpts.useUni = &initSortU;
				initSortU.sortOpts = &sortOpts;
			PreSortMultithreadPipe initSPipe(PIPE_BUFFER_SIZE);
				initSortU.initSPipe = &initSPipe;
			GZipOutStream initOut(0, sortIndOutName.c_str());
				initSortU.initOut = &initOut;
			void* sortThread = startThread(profinmanSortTableThreadFun, &initSortU);
			std::vector<char> tmpStore;
			BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
			BCompTabularReader gfaOut(&blkComp, fastiFN.c_str());
			uintptr_t entryInd = 0;
			while(gfaOut.readNextEntry()){
				tmpStore.clear();
				profinmanSortPackageTableRow(entryInd, &gfaOut, &tmpStore, &fillInds, &fillCode, &fieldSizes);
				initSPipe.writeBytes(&(tmpStore[0]), tmpStore.size());
				entryInd++;
			}
			initSPipe.closeWrite();
			joinThread(sortThread);
		}
	//write to the output
		{
			//open up the output
				std::string ebaseFN(outputName);
				std::string eblockFN = ebaseFN + ".blk";
				std::string efastiFN = ebaseFN + ".tai";
				GZipCompressionMethod ecompMeth;
				BlockCompOutStream eblkComp(0, 0x010000, ebaseFN.c_str(), eblockFN.c_str(), &ecompMeth);
				BCompTabularWriter egfaOut(0, &eblkComp, efastiFN.c_str());
			//open the input
				BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
				BCompTabularReader gfaOut(&blkComp, fastiFN.c_str());
			//open the sorted temp
				std::vector<char> sinItem; sinItem.resize(totItemSize);
				GZipInStream initIn(sortIndOutName.c_str());
			//gogogo
			while(initIn.readBytes(&(sinItem[0]), totItemSize)){
				union {
					uintptr_t asI;
					char asCs[sizeof(uintptr_t)];
				} punmeharder;
				memcpy(punmeharder.asCs, &(sinItem[0]), sizeof(uintptr_t));
				uintptr_t winI = punmeharder.asI;
				gfaOut.readSpecificEntry(winI);
				egfaOut.numEntries = gfaOut.numEntries;
				egfaOut.entrySizes = gfaOut.entrySizes;
				egfaOut.curEntries = gfaOut.curEntries;
				egfaOut.writeNextEntry();
			}
		}
	//clean up
		killFile(sortIndOutName.c_str());
}

ProfinmanSearchSortedTableCells::ProfinmanSearchSortedTableCells(){
	lookTable = 0;
	outputName = 0;
	maxRam = 500000000;
	stdoutName[0] = '-'; stdoutName[1] = 0;
	mySummary = "  Search through a sorted database.";
	myMainDoc = "Usage: profinman findtabs [OPTION] [FILE]*\n"
		"Search through a sorted database.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan findtabs 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".tsv");
		addStringOption("--dump", &outputName, 0, "    Specify the main location to write to.\n    --dump File.tsv\n", &dumpMeta);
	ArgumentParserStrMeta origMeta("Table");
		origMeta.isFile = true;
		origMeta.fileExts.insert(".bctsv");
		addStringOption("--in", &lookTable, 0, "    Specify the table to search.\n    --in File.bctsv\n", &origMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    Specify a target ram usage, in bytes.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserStrVecMeta ocolMeta("Database Column Codes");
		addStringVectorOption("--dcol", &indexCols, 0, "    The column indices to look through in the database.\n    Prefix with an i for integers, an f for floats\n    --dcol i0\n", &ocolMeta);
	ArgumentParserStrVecMeta icolMeta("Search Column Codes");
		addStringVectorOption("--fcol", &lookCols, 0, "    The column indices to look for in the input.\n    Prefix with an i for integers, an f for floats\n    --fcol i0\n", &icolMeta);
}

int ProfinmanSearchSortedTableCells::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	srcFAs.push_back(argv[0]);
	return 1;
}

void ProfinmanSearchSortedTableCells::printExtraGUIInformation(std::ostream* toPrint){
	(*toPrint) << "STRINGVEC\t<|>\tNAME\tInput Files\tFILE\tREAD\t1\t.tsv" << std::endl;
}

ProfinmanSearchSortedTableCells::~ProfinmanSearchSortedTableCells(){
}

int ProfinmanSearchSortedTableCells::posteriorCheck(){
	if(!lookTable || (strlen(lookTable)==0)){
		argumentError = "Need to specify the table to search through.";
		return 1;
	}
	if(!outputName || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "RAM must be positive.";
		return 1;
	}
	const char* cindErr;
	cindErr = profdatatsvCheckColumnSpec(&indexCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCheckColumnSpec(&lookCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCheckConsistent(&indexCols, &lookCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCompileColumnSpec(&indexCols, &fillIndsI, &fillCodeI); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCompileColumnSpec(&lookCols, &fillIndsL, &fillCodeL); if(cindErr){ argumentError = cindErr; return 1; }
	if(srcFAs.size() == 0){
		srcFAs.push_back(stdoutName);
	}
	return 0;
}

#define RECALCULATE_TOTFSIZE totFieldSize = sizeof(uintptr_t); for(uintptr_t ijk = 0; ijk<fieldSizes.size(); ijk++){ totFieldSize += fieldSizes[ijk]; }

#define RESIZE_THINGS \
	resizeCrap.clear();\
	profinmanSortRepackageEntries(&loadedCrap, &resizeCrap, &fillCodeL, &fieldSizes, &changeSizes);\
	memcpy(&(fieldSizes[0]), &(changeSizes[0]), fieldSizes.size()*sizeof(uintptr_t));\
	std::swap(loadedCrap, resizeCrap);\
	RECALCULATE_TOTFSIZE

#define CHECK_NEED_RESIZE(useTSV, useCode, useInd) \
	bool needRS = false;\
	for(uintptr_t i = 0; i<useInd.size(); i++){\
		if(useCode[i] != COLSPEC_STR){ continue; }\
		uintptr_t curEntLen = useTSV.entrySizes[useInd[i]];\
		if(curEntLen > fieldSizes[i]){ changeSizes[i] = curEntLen; needRS = true; }\
	}\
	if(needRS){ RESIZE_THINGS }

void ProfinmanSearchSortedTableCells::runThing(){
	std::vector<uintptr_t> fieldSizes;
	std::vector<uintptr_t> changeSizes;
	uintptr_t totFieldSize;
	for(uintptr_t i = 0; i<fillCodeL.size(); i++){
		if(fillCodeL[i] == COLSPEC_STR){ fieldSizes.push_back(0); }
		else if(fillCodeL[i] == COLSPEC_INT){ fieldSizes.push_back(sizeof(intptr_t)); }
		else{ fieldSizes.push_back(sizeof(double)); }
	}
	RECALCULATE_TOTFSIZE
	changeSizes.insert(changeSizes.end(), fieldSizes.begin(), fieldSizes.end());
	ProfinmanSortTablePiecesUniform sortCompUni;
		sortCompUni.initSPipe = 0;
		sortCompUni.workFolder = 0;
		sortCompUni.initOut = 0;
		sortCompUni.fillCode = &fillCodeL;
		sortCompUni.fieldSizes = &fieldSizes;
	SortOptions sortOpts;
		sortOpts.itemSize = 0;
		sortOpts.maxLoad = maxRam;
		sortOpts.numThread = 1;
		sortOpts.compMeth = profinmanSortTableCompareFun;
		sortOpts.useUni = &sortCompUni;
		sortCompUni.sortOpts = &sortOpts;
	//open the table
	std::string baseFN(lookTable);
	std::string blockFN = baseFN + ".blk";
	std::string fastiFN = baseFN + ".tai";
	GZipCompressionMethod compMeth;
	BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
	BCompTabularReader gfaOut(&blkComp, fastiFN.c_str());
	//run down the files
	std::vector<char> lookCrap;
	std::vector<char> loadedCrap;
	std::vector<char> resizeCrap;
	std::vector<uintptr_t> repackSize;
	std::vector<const char*> repackEnt;
	InStream* curIn = 0;
	OutStream* curOut = 0;
	try{
		curOut = outputName ? (OutStream*)(new FileOutStream(0, outputName)) : (OutStream*)(new ConsoleOutStream()); {
		TSVTabularWriter curOutT(1, curOut);
		uintptr_t item0 = 0;
		uintptr_t curItem = 0;
		for(uintptr_t fi = 0; fi<srcFAs.size(); fi++){
			char* curFN = srcFAs[fi];
			InStream* curIn = ((strcmp(curFN, "-")==0) ? (InStream*)(new ConsoleInStream()) : (InStream*)(new FileInStream(curFN))); {
			TSVTabularReader curInT(1, curIn);
			int hasNextEnt = curInT.readNextEntry();
			while(hasNextEnt || loadedCrap.size()){
				if(hasNextEnt){
					CHECK_NEED_RESIZE(curInT, fillCodeL, fillIndsL)
					//pack into the storage
					profinmanSortPackageTableRow(curItem, &curInT, &loadedCrap, &fillIndsL, &fillCodeL, &fieldSizes);
					curItem++;
					hasNextEnt = curInT.readNextEntry();
				}
				//if enough loaded, do a thing
				if(!hasNextEnt || (loadedCrap.size() >= (uintptr_t)maxRam)){
					//sort the stuff
						sortOpts.itemSize = totFieldSize;
						inMemoryMergesort(curItem - item0, &(loadedCrap[0]), &sortOpts);
					//cocktail binary search
						int curDir = 0;
						uintptr_t dataBLowInd = 0;
						uintptr_t dataBHigInd = gfaOut.getNumEntries();
						uintptr_t lowSearchInd = 0;
						uintptr_t higSearchInd = loadedCrap.size() / totFieldSize;
						while(lowSearchInd < higSearchInd){
							uintptr_t curLookInd = curDir ? (higSearchInd-1) : lowSearchInd;
							//prepare the output vectors
							char strIndex[4*sizeof(uintptr_t)+4];
							char* curLookPack = &(loadedCrap[curLookInd*totFieldSize]);
							{
								union {
									uintptr_t asI;
									char asCs[sizeof(uintptr_t)];
								} punmeharder;
								memcpy(punmeharder.asCs, curLookPack, sizeof(uintptr_t));
								uintptr_t repIndex = punmeharder.asI;
								sprintf(strIndex, "%ju", (uintmax_t)repIndex);
							}
							repackSize.clear(); repackSize.push_back(strlen(strIndex));
							repackEnt.clear(); repackEnt.push_back(strIndex);
							//lower bound search
							uintptr_t lowBLowI = dataBLowInd;
							{
								uintptr_t count = dataBHigInd - lowBLowI;
								while(count){
									//get the entry
									uintptr_t step = count / 2;
									uintptr_t it = lowBLowI + step;
									gfaOut.readSpecificEntry(it);
									CHECK_NEED_RESIZE(gfaOut, fillCodeI, fillIndsI)
									sortOpts.itemSize = totFieldSize;
									lookCrap.clear();
									profinmanSortPackageTableRow(0, &gfaOut, &lookCrap, &fillIndsI, &fillCodeI, &fieldSizes);
									//do the comparison
									if(sortOpts.compMeth(sortOpts.useUni, &(lookCrap[0]), curLookPack)){
										lowBLowI = it + 1;
										count = count - (step + 1);
									}
									else{
										count = step;
									}
								}
							}
							//upper bound search
							uintptr_t higBLowI = lowBLowI;
							{
								uintptr_t count = dataBHigInd - higBLowI;
								while(count){
									//get the entry
									uintptr_t step = count / 2;
									uintptr_t it = higBLowI + step;
									gfaOut.readSpecificEntry(it);
									CHECK_NEED_RESIZE(gfaOut, fillCodeI, fillIndsI)
									sortOpts.itemSize = totFieldSize;
									lookCrap.clear();
									profinmanSortPackageTableRow(0, &gfaOut, &lookCrap, &fillIndsI, &fillCodeI, &fieldSizes);
									//do the comparison
									if(!sortOpts.compMeth(sortOpts.useUni, curLookPack, &(lookCrap[0]))){
										higBLowI = it + 1;
										count = count - (step + 1);
									}
									else{
										count = step;
									}
								}
							}
							//report
							for(uintptr_t ri = lowBLowI; ri < higBLowI; ri++){
								gfaOut.readSpecificEntry(ri);
								repackSize.resize(1); repackSize.insert(repackSize.end(), gfaOut.entrySizes, gfaOut.entrySizes + gfaOut.numEntries);
								repackEnt.resize(1); repackEnt.insert(repackEnt.end(), gfaOut.curEntries, gfaOut.curEntries + gfaOut.numEntries);
								curOutT.numEntries = repackSize.size();
								curOutT.entrySizes = &(repackSize[0]);
								curOutT.curEntries = &(repackEnt[0]);
								curOutT.writeNextEntry();
							}
							//update direction
							if(curDir){ dataBHigInd = higBLowI; higSearchInd--; }
								else{ dataBLowInd = lowBLowI; lowSearchInd++; }
							curDir = !curDir;
						}
					//prepare for the next round
						loadedCrap.clear();
						item0 = curItem;
				}
			}
			} delete(curIn); curIn = 0;
		}
		} delete(curOut); curOut = 0;
	}
	catch(std::exception& err){
		if(curIn){ delete(curIn); }
		if(curOut){ delete(curOut); }
		throw;
	}
}

ProfinmanSearchTableCells::ProfinmanSearchTableCells(){
	lookTable = 0;
	outputName = 0;
	maxRam = 500000000;
	stdoutName[0] = '-'; stdoutName[1] = 0;
	mySummary = "  Search through a database.";
	myMainDoc = "Usage: profinman findtab [OPTION] [FILE]*\n"
		"Search through a database.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan findtab 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".tsv");
		addStringOption("--dump", &outputName, 0, "    Specify the main location to write to.\n    --dump File.tsv\n", &dumpMeta);
	ArgumentParserStrMeta origMeta("Table");
		origMeta.isFile = true;
		origMeta.fileExts.insert(".bctsv");
		addStringOption("--in", &lookTable, 0, "    Specify the table to search.\n    --in File.bctsv\n", &origMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    Specify a target ram usage, in bytes.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserStrVecMeta ocolMeta("Database Column Codes");
		addStringVectorOption("--dcol", &indexCols, 0, "    The column indices to look through in the database.\n    Prefix with an i for integers, an f for floats\n    --dcol i0\n", &ocolMeta);
	ArgumentParserStrVecMeta icolMeta("Search Column Codes");
		addStringVectorOption("--fcol", &lookCols, 0, "    The column indices to look for in the input.\n    Prefix with an i for integers, an f for floats\n    --fcol i0\n", &icolMeta);
}

int ProfinmanSearchTableCells::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	srcFAs.push_back(argv[0]);
	return 1;
}

void ProfinmanSearchTableCells::printExtraGUIInformation(std::ostream* toPrint){
	(*toPrint) << "STRINGVEC\t<|>\tNAME\tInput Files\tFILE\tREAD\t1\t.tsv" << std::endl;
}

ProfinmanSearchTableCells::~ProfinmanSearchTableCells(){
}

int ProfinmanSearchTableCells::posteriorCheck(){
	if(!lookTable || (strlen(lookTable)==0)){
		argumentError = "Need to specify the table to search through.";
		return 1;
	}
	if(!outputName || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "RAM must be positive.";
		return 1;
	}
	const char* cindErr;
	cindErr = profdatatsvCheckColumnSpec(&indexCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCheckColumnSpec(&lookCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCheckConsistent(&indexCols, &lookCols); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCompileColumnSpec(&indexCols, &fillIndsI, &fillCodeI); if(cindErr){ argumentError = cindErr; return 1; }
	cindErr = profdatatsvCompileColumnSpec(&lookCols, &fillIndsL, &fillCodeL); if(cindErr){ argumentError = cindErr; return 1; }
	if(srcFAs.size() == 0){
		srcFAs.push_back(stdoutName);
	}
	return 0;
}

void ProfinmanSearchTableCells::runThing(){
	std::vector<uintptr_t> fieldSizes;
	std::vector<uintptr_t> changeSizes;
	uintptr_t totFieldSize;
	for(uintptr_t i = 0; i<fillCodeL.size(); i++){
		if(fillCodeL[i] == COLSPEC_STR){ fieldSizes.push_back(0); }
		else if(fillCodeL[i] == COLSPEC_INT){ fieldSizes.push_back(sizeof(intptr_t)); }
		else{ fieldSizes.push_back(sizeof(double)); }
	}
	RECALCULATE_TOTFSIZE
	changeSizes.insert(changeSizes.end(), fieldSizes.begin(), fieldSizes.end());
	ProfinmanSortTablePiecesUniform sortCompUni;
		sortCompUni.initSPipe = 0;
		sortCompUni.workFolder = 0;
		sortCompUni.initOut = 0;
		sortCompUni.fillCode = &fillCodeL;
		sortCompUni.fieldSizes = &fieldSizes;
	SortOptions sortOpts;
		sortOpts.itemSize = 0;
		sortOpts.maxLoad = maxRam;
		sortOpts.numThread = 1;
		sortOpts.compMeth = profinmanSortTableCompareFun;
		sortOpts.useUni = &sortCompUni;
		sortCompUni.sortOpts = &sortOpts;
	std::string baseFN(lookTable);
	std::string blockFN = baseFN + ".blk";
	std::string fastiFN = baseFN + ".tai";
	GZipCompressionMethod compMeth;
	//run down the files
	std::vector<char> lookCrap;
	std::vector<char> loadedCrap;
	std::vector<char> resizeCrap;
	std::vector<uintptr_t> repackSize;
	std::vector<const char*> repackEnt;
	InStream* curIn = 0;
	OutStream* curOut = 0;
	try{
		curOut = outputName ? (OutStream*)(new FileOutStream(0, outputName)) : (OutStream*)(new ConsoleOutStream()); {
		TSVTabularWriter curOutT(1, curOut);
		uintptr_t item0 = 0;
		uintptr_t curItem = 0;
		for(uintptr_t fi = 0; fi<srcFAs.size(); fi++){
			char* curFN = srcFAs[fi];
			InStream* curIn = ((strcmp(curFN, "-")==0) ? (InStream*)(new ConsoleInStream()) : (InStream*)(new FileInStream(curFN))); {
			TSVTabularReader curInT(1, curIn);
			int hasNextEnt = curInT.readNextEntry();
			while(hasNextEnt || loadedCrap.size()){
				if(hasNextEnt){
					CHECK_NEED_RESIZE(curInT, fillCodeL, fillIndsL)
					//pack into the storage
					profinmanSortPackageTableRow(curItem, &curInT, &loadedCrap, &fillIndsL, &fillCodeL, &fieldSizes);
					curItem++;
					hasNextEnt = curInT.readNextEntry();
				}
				//if enough loaded, do a thing
				if(!hasNextEnt || (loadedCrap.size() >= (uintptr_t)maxRam)){
					//sort the stuff
						sortOpts.itemSize = totFieldSize;
						inMemoryMergesort(curItem - item0, &(loadedCrap[0]), &sortOpts);
					//open the table
						BlockCompInStream blkComp(baseFN.c_str(), blockFN.c_str(), &compMeth);
						BCompTabularReader gfaOut(&blkComp, fastiFN.c_str());
					//go entry by entry
						while(gfaOut.readNextEntry()){
							CHECK_NEED_RESIZE(gfaOut, fillCodeI, fillIndsI)
							sortOpts.itemSize = totFieldSize;
							lookCrap.clear();
							profinmanSortPackageTableRow(0, &gfaOut, &lookCrap, &fillIndsI, &fillCodeI, &fieldSizes);
							//package the thing for a join
							repackSize.clear(); repackSize.push_back(0); repackSize.insert(repackSize.end(), gfaOut.entrySizes, gfaOut.entrySizes + gfaOut.numEntries);
							repackEnt.clear(); repackEnt.push_back(0); repackEnt.insert(repackEnt.end(), gfaOut.curEntries, gfaOut.curEntries + gfaOut.numEntries);
							//figure the bounds of which loaded items could match
							char* lowMat = whodunSortLowerBound(loadedCrap.size() / totFieldSize, &(loadedCrap[0]), &(lookCrap[0]), &sortOpts);
							char* higMat = whodunSortUpperBound(loadedCrap.size() / totFieldSize, &(loadedCrap[0]), &(lookCrap[0]), &sortOpts);
							//see if any of them match
							for(char* curLook = lowMat; curLook < higMat; curLook += totFieldSize){
								uintptr_t curIndex = item0 + ((curLook - &(loadedCrap[0]))/totFieldSize);
								char strIndex[4*sizeof(uintptr_t)+4];
								sprintf(strIndex, "%ju", (uintmax_t)curIndex);
								repackSize[0] = strlen(strIndex);
								repackEnt[0] = strIndex;
								curOutT.numEntries = repackSize.size();
								curOutT.entrySizes = &(repackSize[0]);
								curOutT.curEntries = &(repackEnt[0]);
								curOutT.writeNextEntry();
							}
						}
					//prepare for the next round
						loadedCrap.clear();
						item0 = curItem;
				}
			}
			} delete(curIn); curIn = 0;
		}
		} delete(curOut); curOut = 0;
	}
	catch(std::exception& err){
		if(curIn){ delete(curIn); }
		if(curOut){ delete(curOut); }
		throw;
	}
}

ProfinmanSortTableSearchResult::ProfinmanSortTableSearchResult(){
	lookTable = 0;
	maxRam = 500000000;
	numThread = 1;
	workFolder = 0;
	mySummary = "  Sort database search results by initial record.";
	myMainDoc = "Usage: profinman sortfindtab [OPTION] [FILE]*\n"
		"Sort databasesearch results by initial record.\n"
		"The OPTIONS are:\n";
	myVersionDoc = "ProFinMan sortfindtab 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta dumpMeta("Main Output Location");
		dumpMeta.isFile = true;
		dumpMeta.fileWrite = true;
		dumpMeta.fileExts.insert(".tsv");
		addStringOption("--dump", &outputName, 0, "    Specify the main location to write to.\n    --dump File.tsv\n", &dumpMeta);
	ArgumentParserStrMeta origMeta("Unsorted");
		origMeta.isFile = true;
		origMeta.fileExts.insert(".tsv");
		addStringOption("--in", &lookTable, 0, "    Specify the results to sort.\n    --in File.tsv\n", &origMeta);
	ArgumentParserIntMeta ramMeta("RAM Usage");
		addIntegerOption("--ram", &maxRam, 0, "    Specify a target ram usage, in bytes.\n    --ram 500000000\n", &ramMeta);
	ArgumentParserIntMeta threadMeta("Threads");
		addIntegerOption("--thread", &numThread, 0, "    How many threads to use.\n    --thread 1\n", &threadMeta);
	ArgumentParserStrMeta workMeta("Working Folder");
		workMeta.isFolder = true;
		addStringOption("--work", &workFolder, 0, "    The folder to put temporary files in.\n    --work Folder\n", &workMeta);
}

ProfinmanSortTableSearchResult::~ProfinmanSortTableSearchResult(){
}

int ProfinmanSortTableSearchResult::posteriorCheck(){
	if(!lookTable || (strlen(lookTable)==0)){
		lookTable = 0;
	}
	if(!outputName || (strlen(outputName)==0)){
		outputName = 0;
	}
	if(maxRam <= 0){
		argumentError = "RAM must be positive.";
		return 1;
	}
	if(numThread <= 0){
		argumentError = "Thread count must be positive.";
		return 1;
	}
	if((workFolder == 0) || (strlen(workFolder)==0)){
		argumentError = "Need to specify a working folder.";
		return 1;
	}
	return 0;
}

void ProfinmanSortTableSearchResult::runThing(){
	GZipCompressionMethod compMeth;
	//names of temporary files
		std::string bctabOutName = workFolder;
			bctabOutName.append(pathElementSep);
			bctabOutName.append("dsres_tmp.bctsv");
		std::string bctabBOutName = bctabOutName + ".blk";
		std::string bctabIOutName = bctabOutName + ".tai";
		std::string sbctabOutName = workFolder;
			sbctabOutName.append(pathElementSep);
			sbctabOutName.append("dsres_sort_tmp.bctsv");
		std::string sbctabBOutName = sbctabOutName + ".blk";
		std::string sbctabIOutName = sbctabOutName + ".tai";
	//drain the thing to the work folder
		{
			BlockCompOutStream eblkComp(0, 0x010000, bctabOutName.c_str(), bctabBOutName.c_str(), &compMeth);
			BCompTabularWriter egfaOut(0, &eblkComp, bctabIOutName.c_str());
			InStream* txtTabIn = 0;
			if(lookTable){ txtTabIn = new FileInStream(lookTable); } else{ txtTabIn = new ConsoleInStream(); }
			try{
				TSVTabularReader txtTabT(0, txtTabIn);
				while(txtTabT.readNextEntry()){
					egfaOut.numEntries = txtTabT.numEntries;
					egfaOut.entrySizes = txtTabT.entrySizes;
					egfaOut.curEntries = txtTabT.curEntries;
					egfaOut.writeNextEntry();
				}
			}
			catch(...){
				if(txtTabIn){ delete(txtTabIn); }
				throw;
			}
			delete(txtTabIn);
		}
	//use the other sorting program to actually sort
		std::vector<char> lookTTmp(bctabOutName.begin(), bctabOutName.end());
		std::vector<char> dumpTTmp(sbctabOutName.begin(), sbctabOutName.end());
		char indSpecTmp[3]; indSpecTmp[0] = 'i'; indSpecTmp[1] = '0'; indSpecTmp[2] = 0;
		ProfinmanSortTableCells subTask;
			subTask.lookTable = &(lookTTmp[0]);
			subTask.outputName = &(dumpTTmp[0]);
			subTask.indexCols.push_back(indSpecTmp);
			subTask.maxRam = maxRam;
			subTask.numThread = numThread;
			subTask.workFolder = workFolder;
		if(subTask.posteriorCheck()){
			throw std::runtime_error("...what?");
		}
		subTask.runThing();
	//dump table
		{
			BlockCompInStream blkComp(sbctabOutName.c_str(), sbctabBOutName.c_str(), &compMeth);
			BCompTabularReader gfaOut(&blkComp, sbctabIOutName.c_str());
			OutStream* txtTabOut = 0;
			if(outputName){ txtTabOut = new FileOutStream(0, outputName); }else{ txtTabOut = new ConsoleOutStream(); }
			try{
				TSVTabularWriter txtTabT(0, txtTabOut);
				while(gfaOut.readNextEntry()){
					txtTabT.numEntries = gfaOut.numEntries;
					txtTabT.entrySizes = gfaOut.entrySizes;
					txtTabT.curEntries = gfaOut.curEntries;
					txtTabT.writeNextEntry();
				}
			}
			catch(...){
				if(txtTabOut){ delete(txtTabOut); }
				throw;
			}
			delete(txtTabOut);
		}
	//clean up
		killFile(bctabOutName.c_str());
		killFile(bctabBOutName.c_str());
		killFile(bctabIOutName.c_str());
		killFile(sbctabOutName.c_str());
		killFile(sbctabBOutName.c_str());
		killFile(sbctabIOutName.c_str());
}



