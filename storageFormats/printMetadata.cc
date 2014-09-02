
#include <algorithm>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "reportingUtilsEnvironment.h"
#include "reportingUtils.hpp"
#include "eventFileWriter.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "print metadata information for a rootpwa file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -i inputFile [-r]" << endl
	     << "    where:" << endl
	     << "        -i file    input file in ROOTPWA format (default: ./input.root)" << endl
	     << "        -r         recalculate hash and compare it with the stored one (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}

int main(int argc, char** argv)
{

	printCompilerInfo();
	printLibraryInfo();
	printGitHash();
	cout << endl;

	string dataMetadataName = "dataMetadata";
	string inTreeName = "rootPwaEvtTree";
	string prodKinMomentaLeafName = "prodKinMomenta";
	string decayKinMomentaLeafName = "decayKinMomenta";

	const string progName = argv[0];
	string inputFileName  = "input.root";
	bool recalculateHash = false;

	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");

	extern char* optarg;
	int c;
	while((c = getopt(argc, argv, "i:rh")) != -1)
	{
		switch(c) {
		case 'i':
			inputFileName = optarg;
			break;
		case 'r':
			recalculateHash = true;
			break;
		case 'h':
			usage(progName);
			break;
		}
	}
	TFile* inputFile = TFile::Open(inputFileName.c_str(), "READ");
	if(not inputFile) {
		printErr << "could not open input file. Aborting..." << endl;
		return 1;
	}
	eventMetadata* metadata = (eventMetadata*)inputFile->Get(dataMetadataName.c_str());
	if(metadata) {
		// we are reading a datafile
		TTree* inputTree = (TTree*)inputFile->Get(inTreeName.c_str());
		if(not inputTree) {
			printErr << "reading a datafile but did not find event tree "
			         << "with name '" << inTreeName << "'. Aborting..." << endl;
			return 1;
		}
		const vector<string>& additionalVariableNames = metadata->additionalSavedVariableLables();
		{
			TObjArray* branchList = inputTree->GetListOfBranches();
			for(int i = 0; i < branchList->GetEntries(); ++i) {
				TBranch* branch = (TBranch*)(*branchList)[i];
				const string branchName = branch->GetName();
				if((branchName == prodKinMomentaLeafName) or (branchName == decayKinMomentaLeafName))
				{
					continue;
				}
				if(find(additionalVariableNames.begin(), additionalVariableNames.end(), branchName) == additionalVariableNames.end()) {
					printErr << "additional branch '" << branchName << "' present in metadata, but not found in tree. Aborting..." << endl;
					return 1;
				}
			}
		}
		if(recalculateHash) {
			printInfo << "recalculating hash..." << endl;
			const string calculatedHash = eventFileWriter::calculateHash(inputTree, additionalVariableNames, true);
			if(calculatedHash != metadata->contentHash()) {
				printErr << "hash verification failed, hash from metadata '" << metadata->contentHash() << "' does "
				         << "not match with calculated hash '" << calculatedHash << "'. Aborting..." << endl;
			} else {
				cout << endl;
				printSucc << "recalculated hash matches with hash from metadata." << endl;
				cout << endl;
			}
		}
		printInfo << *metadata << endl;
		printInfo << "additional information:" << endl
		          << "    number of events in file ... " << inputTree->GetEntries() << endl
		          << "    additional branches ........ " << additionalVariableNames << endl;
		cout << endl;
		return 0;
	}
	// possibilities to read other file formats, i.e. amplitudes, fit inputs, fit outputs, etc.
	printErr << "could not find any metadata." << endl;
	return 1;
}
