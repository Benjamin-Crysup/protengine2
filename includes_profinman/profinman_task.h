#ifndef PROFINMAN_TASK_H
#define PROFINMAN_TASK_H 1

#include <vector>
#include <fstream>
#include <utility>
#include <stdio.h>
#include <stdint.h>

#include "whodun_args.h"

/**A profinman action.*/
class ProfinmanAction : public ArgumentParser{
public:
	/**Let subclasses handle.*/
	virtual ~ProfinmanAction();
	/**
	 * Run the action.
	 */
	virtual void runThing() = 0;
	/**The summary of this action type.*/
	const char* mySummary;
};

//************************************************************************
//SEQUENCE STUFF
//************************************************************************

/**Use to prepare a sequence for use.*/
class ProfinmanBlockSequence : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanBlockSequence();
	~ProfinmanBlockSequence();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void runThing();
	
	/**The base name of the dumps.*/
	char* dumpBaseName;
	/**The files to block.*/
	std::vector<const char*> srcFAs;
};

/**Dump sequence to fasta*/
class ProfinmanDumpSequence : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanDumpSequence();
	~ProfinmanDumpSequence();
	int posteriorCheck();
	void runThing();
	
	/**The base name of the thing to dump.*/
	char* dumpBaseName;
	/**The name of the file to write to.*/
	char* outputName;
};

/**Direct search a sequence.*/
class ProfinmanSearchSequence : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSearchSequence();
	/**The base name of the reference.*/
	char* dumpBaseName;
	/**The name of the search fasta.*/
	char* searchName;
	/**The name of the file to write to.*/
	char* outputName;
	/**The maximum number of bytes of search sequence to load.*/
	intptr_t maxRam;
	/**Output results in text.*/
	bool txtOut;
	int posteriorCheck();
	void runThing();
};

//************************************************************************
//MATCH EXAMINATION
//************************************************************************

/**Get sequence near a match.*/
class ProfinmanGetMatchRegion : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanGetMatchRegion();
	/**The base name of the reference.*/
	char* dumpBaseName;
	/**The name of the search result file.*/
	char* matchName;
	/**The number of upstream bases to get.*/
	intptr_t numPre;
	/**The number of downstream bases to get.*/
	intptr_t numPost;
	/**The place to write the output.*/
	char* outputName;
	int posteriorCheck();
	void runThing();
};

/**Sort binary search results*/
class ProfinmanSortSearchResults : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSortSearchResults();
	~ProfinmanSortSearchResults();
	int posteriorCheck();
	void runThing();
	/**The base name of the thing to dump.*/
	char* origSRName;
	/**The maximum number of bytes to use.*/
	intptr_t maxRam;
	/**The number of threads to use.*/
	intptr_t numThread;
	/**The place to write the output.*/
	char* outputName;
	/**The folder to work in.*/
	char* workFolder;
};

/**Get the name of the sequence that was matched.*/
class ProfinmanGetMatchName : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanGetMatchName();
	int posteriorCheck();
	void runThing();
	/**The base name of the reference.*/
	char* dumpBaseName;
	/**The name of the search result file.*/
	char* matchName;
	/**The place to write the output.*/
	char* outputName;
	/**The maximum number of bytes to use.*/
	intptr_t maxRam;
};

/**Filter hits that are valid (according to a digest file).*/
class ProfinmanFilterDigestMatches : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanFilterDigestMatches();
	int posteriorCheck();
	void runThing();
	/**The base name of the reference.*/
	char* dumpBaseName;
	/**The name of the search result file.*/
	char* matchName;
	/**The digest file to compare to.*/
	char* digestName;
	/**The place to write the output.*/
	char* outputName;
	/**The maximum number of bytes to use.*/
	intptr_t maxRam;
};

//TODO

//************************************************************************
//SUFFIX ARRAYS
//************************************************************************

/**Use to build a suffix array.*/
class ProfinmanBuildReference : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanBuildReference();
	/**Clean up.*/
	~ProfinmanBuildReference();
	/**The base name of the reference.*/
	char* referenceName;
	/**The base name of the combo file.*/
	char* comboName;
	/**The maximum number of bytes to use.*/
	intptr_t maxRam;
	/**The number of threads to use.*/
	intptr_t numThread;
	/**The folder to work in.*/
	char* workFolder;
	/**The file to recover with.*/
	char* recoverFile;
	
	int posteriorCheck();
	void runThing();
	
	/**The opened recovery file*/
	std::ofstream* recoverStream = 0;
};

/**Use to dump a suffix array.*/
class ProfinmanDumpReference : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanDumpReference();
	/**Clean up.*/
	~ProfinmanDumpReference();
	/**The base name of the combo file.*/
	char* comboName;
	/**The base name of the reference.*/
	char* referenceName;
	/**The number of downstream bases to get.*/
	intptr_t numPost;
	int posteriorCheck();
	void runThing();
};

/**Direct search a suffix array.*/
class ProfinmanSearchReference : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSearchReference();
	/**The base name of the reference.*/
	char* referenceName;
	/**The base name of the combo file.*/
	char* comboName;
	/**The name of the file to write to.*/
	char* outputName;
	/**The name of the search fasta.*/
	char* searchName;
	/**Output results in text.*/
	bool txtOut;
	int posteriorCheck();
	void runThing();
};

/**Use to merge suffix arrays.*/
class ProfinmanMergeReference : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanMergeReference();
	/**Clean up.*/
	~ProfinmanMergeReference();
	/**The base name of the first reference.*/
	char* refAName;
	/**The base name of the second reference.*/
	char* refBName;
	/**The base name of the first combo.*/
	char* comAName;
	/**The base name of the second combo.*/
	char* comBName;
	/**The base name of the output reference.*/
	char* referenceName;
	/**The base name of the output combo file.*/
	char* comboName;
	
	int posteriorCheck();
	void runThing();
};

//TODO

//************************************************************************
//DATABASES
//************************************************************************

/**Pack a tsv into a database.*/
class ProfinmanPackTable : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanPackTable();
	~ProfinmanPackTable();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void runThing();
	/**The base name of the thing to dump.*/
	char* dumpBaseName;
	/**The files to block.*/
	std::vector<const char*> srcFAs;
	
	/**Reliable storage for the name for stdout/stdin*/
	char stdoutName[2];
};

/**Sort a table.*/
class ProfinmanSortTableCells : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSortTableCells();
	~ProfinmanSortTableCells();
	int posteriorCheck();
	void runThing();
	/**The table to look through.*/
	char* lookTable;
	/**The file to write to.*/
	char* outputName;
	/**The database columns containing the IDs.*/
	std::vector<char*> indexCols;
	/**The maximum number of bytes of search sequence to load.*/
	intptr_t maxRam;
	/**The number of threads to use.*/
	intptr_t numThread;
	/**The folder to work in.*/
	char* workFolder;
	
	/**Reliable storage for the name for stdout/stdin*/
	char stdoutName[2];
	/**Store compiled indices.*/
	std::vector<intptr_t> fillInds;
	/**Store compiled indices.*/
	std::vector<int> fillCode;
};

/**Search through a sorted binary table.*/
class ProfinmanSearchSortedTableCells : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSearchSortedTableCells();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	~ProfinmanSearchSortedTableCells();
	int posteriorCheck();
	void runThing();
	/**The table to look through.*/
	char* lookTable;
	/**The files to find for.*/
	std::vector<char*> srcFAs;
	/**The file to write to.*/
	char* outputName;
	/**The input columns containing the IDs.*/
	std::vector<char*> lookCols;
	/**The database columns containing the IDs.*/
	std::vector<char*> indexCols;
	/**The maximum number of bytes of search sequence to load.*/
	intptr_t maxRam;
	
	/**Reliable storage for the name for stdout/stdin*/
	char stdoutName[2];
	/**Store compiled indices.*/
	std::vector<intptr_t> fillIndsL;
	/**Store compiled indices.*/
	std::vector<int> fillCodeL;
	/**Store compiled indices.*/
	std::vector<intptr_t> fillIndsI;
	/**Store compiled indices.*/
	std::vector<int> fillCodeI;
};

/**Search through an unsorted binary table.*/
class ProfinmanSearchTableCells : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSearchTableCells();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	~ProfinmanSearchTableCells();
	int posteriorCheck();
	void runThing();
	/**The table to look through.*/
	char* lookTable;
	/**The files to find for.*/
	std::vector<char*> srcFAs;
	/**The file to write to.*/
	char* outputName;
	/**The input columns containing the IDs.*/
	std::vector<char*> lookCols;
	/**The database columns containing the IDs.*/
	std::vector<char*> indexCols;
	/**The maximum number of bytes of search sequence to load.*/
	intptr_t maxRam;
	
	/**Reliable storage for the name for stdout/stdin*/
	char stdoutName[2];
	/**Store compiled indices.*/
	std::vector<intptr_t> fillIndsL;
	/**Store compiled indices.*/
	std::vector<int> fillCodeL;
	/**Store compiled indices.*/
	std::vector<intptr_t> fillIndsI;
	/**Store compiled indices.*/
	std::vector<int> fillCodeI;
};

/**Sort the joined results of a table search by search record.*/
class ProfinmanSortTableSearchResult : public ProfinmanAction{
public:
	/**Set up an empty action.*/
	ProfinmanSortTableSearchResult();
	~ProfinmanSortTableSearchResult();
	int posteriorCheck();
	void runThing();
	/**The results to sort.*/
	char* lookTable;
	/**The file to write to.*/
	char* outputName;
	/**The maximum number of bytes of search sequence to load.*/
	intptr_t maxRam;
	/**The number of threads to use.*/
	intptr_t numThread;
	/**The folder to work in.*/
	char* workFolder;
};

//TODO

//************************************************************************
//RANDOM CRAP
//************************************************************************

/**Essentially memcmp.*/
bool memBlockCompare(const std::pair<const char*,uintptr_t>& itemA, const std::pair<const char*,uintptr_t>& itemB);

/**Do a comparison where it is known that all characters up to a certain point are the same.*/
class SingleCharMemBlockCompare{
public:
	/**The index to compare at.*/
	uintptr_t compInd;
	bool operator ()(const std::pair<const char*,uintptr_t>& itemA, const std::pair<const char*,uintptr_t>& itemB);
};

/**The size of a search result entry (4 numbers)*/
#define MATCH_ENTRY_SIZE 32

/**
 * Output a search result to a file.
 * @param lookFor The sequence index that was found.
 * @param foundIn The reference sequence it was found in.
 * @param foundAt THe location in the reference it was found at.
 * @param foundTo The location after the end in the reference.
 * @param outputTo The file to write to.
 * @param asText Write as text (or binary).
 */
void outputSearchResult(uintptr_t lookFor, uintptr_t foundIn, uintptr_t foundAt, uintptr_t foundTo, FILE* outputTo, bool asText);

/**The size of an entry in a working combo file.*/
#define COMBO_SORT_ENTRY_SIZE 32
/**The size of an entry in a finalized combo file.*/
#define COMBO_ENTRY_SIZE 16

#endif
