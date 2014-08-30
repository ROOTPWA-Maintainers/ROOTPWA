
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "reportingUtilsEnvironment.h"
#include "reportingUtils.hpp"
#include "dataFileWriter.h"


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
	     << " -i inputFile" << endl
	     << "    where:" << endl
	     << "        -i file    input file in ROOTPWA format (default: ./input.root)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}

int main(int argc, char** argv)
{
	// TODO: Add option to recalculate and thereby verify the contentHash

	printCompilerInfo();
	printLibraryInfo();
	printGitHash();
	cout << endl;

	string metadataName = "dataMetadata";
	string inTreeName = "rootPwaEvtTree";
	string prodKinMomentaLeafName = "prodKinMomenta";
	string decayKinMomentaLeafName = "decayKinMomenta";

	const string progName = argv[0];
	string inputFileName  = "input.root";

	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");

	extern char* optarg;
	int c;
	while((c = getopt(argc, argv, "i:h")) != -1)
	{
		switch(c) {
		case 'i':
			inputFileName = optarg;
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
	dataMetadata* metadata = (dataMetadata*)inputFile->Get(metadataName.c_str());
	if(metadata) {
		// we are reading a datafile
		TTree* inputTree = (TTree*)inputFile->Get(inTreeName.c_str());
		if(not inputTree) {
			printErr << "reading a datafile but did not find event tree "
			         << "with name '" << inTreeName << "'. Aborting..." << endl;
			return 1;
		}
		std::vector<string> additionalBranches;
		TObjArray* branchList = inputTree->GetListOfBranches();
		for(int i = 0; i < branchList->GetEntries(); ++i) {
			TBranch* branch = (TBranch*)(*branchList)[i];
			string branchName = branch->GetName();
			if((branchName == prodKinMomentaLeafName) or (branchName == decayKinMomentaLeafName))
			{
				continue;
			}
			additionalBranches.push_back(branch->GetName());
		}
		printInfo << *metadata << endl;
		printInfo << "additional information:" << endl
		          << "    number of events in file ... " << inputTree->GetEntries() << endl
		          << "    additional branches ........ " << additionalBranches << endl;
		cout << endl;
		return 0;
	}
	// possibilities to read other file formats, i.e. amplitudes, fit inputs, fit outputs, etc.
	printErr << "could not find any metadata." << endl;
	return 1;
}
