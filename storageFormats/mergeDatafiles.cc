
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

	cerr << "merge several datafiles into one" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-f] outputFile inputFile1 inputFile2 ..." << endl
	     << "    where:" << endl
	     << "        -f         overwrite output file if it exists" << endl
	     << endl;
	exit(errCode);
}

int main(int argc, char** argv)
{

	printCompilerInfo();
	printLibraryInfo();
	printGitHash();
	cout << endl;

	const string progName = argv[0];
	bool force = false;

	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");

//	extern char* optarg;
	int c;
	while((c = getopt(argc, argv, "fh")) != -1)
	{
		switch(c) {
		case 'f':
			force = true;
			break;
		case 'h':
			usage(progName);
			break;
		}
	}
	string outputFileName = "";
	vector<string> inputFileNames;
	{
		unsigned int startingPoint = 1;
		if(force) {
			++startingPoint;
		}
		outputFileName = argv[startingPoint++];
		for(int i = startingPoint; i < argc; ++i) {
			inputFileNames.push_back(argv[i]);
		}
	}
	printInfo << "target file: " << outputFileName << endl;
	for(unsigned int i = 0; i < inputFileNames.size(); ++i) {
		printInfo << "source file " << i << ": " << inputFileNames[i] << endl;
	}
	string outputFileOpenOption = "NEW";
	if(force) {
		outputFileOpenOption = "RECREATE";
	}
	TFile* outputFile = TFile::Open(outputFileName.c_str(), outputFileOpenOption.c_str());
	if(not outputFile) {
		printErr << "could not open output file '" << outputFileName << "'. Aborting..." << endl;
		return 1;
	}
	vector<TFile*> inputFiles(inputFileNames.size(), 0);
	for(unsigned int i = 0; i < inputFileNames.size(); ++i) {
		inputFiles[i] = TFile::Open(inputFileNames[i].c_str(), "READ");
		if(not inputFiles[i]) {
			printErr << "could not open input file '" << inputFileNames[i] << "'. Aborting..." << endl;
			return 1;
		}
	}
	vector<pair<const eventMetadata*, TTree*> > inputData(inputFiles.size(), pair<const eventMetadata*, TTree*>(0, 0));
	outputFile->cd();
	eventMetadata metadata;
	for(unsigned int i = 0; i < inputFiles.size(); ++i) {
		eventMetadata* inputMetadata = (eventMetadata*)inputFiles[i]->Get(eventMetadata::objectNameInFile.c_str());
		TTree* inputTree = (TTree*)inputFiles[i]->Get(eventMetadata::eventTreeName.c_str());
		inputData[i] = pair<const eventMetadata*, TTree* >(inputMetadata, inputTree);
	}
	TTree* outputTree = metadata.merge(inputData);
	if(not outputTree) {
		printErr << "merge failed. Aborting..." << endl;
		return 1;
	}
	outputFile->cd();
	metadata.Write(eventMetadata::objectNameInFile.c_str());
	outputTree->Write();
	outputFile->Close();
	for(unsigned int i = 0; i < inputFiles.size(); ++i) {
		inputFiles[i]->Close();
	}
	return 0;
}
