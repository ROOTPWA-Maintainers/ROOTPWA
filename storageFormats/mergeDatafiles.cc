
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
	     << " [-a -f] outputFile inputFile1 inputFile2 ..." << endl
	     << "    where:" << endl
	     << "        -a         accept different metadata and merge to combined bin" << endl
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
	bool mergeDiffMeta = false;
	bool force         = false;

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");
#endif

//	extern char* optarg;
	int c;
	extern int optind;
	while((c = getopt(argc, argv, "afh")) != -1)
	{
		switch(c) {
		case 'a':
			mergeDiffMeta = true;
			break;
		case 'f':
			force = true;
			break;
		case 'h':
			usage(progName);
			break;
		}
	}
	if (argc - optind < 3) {
		printErr << "you have to specify at least two data files to be merged. Aborting..." << endl;
		usage(progName, 1);
	}

	vector<string> inputFileNames;
	string outputFileName = argv[optind];
	for(int i = optind+1; i < argc; ++i) {
		inputFileNames.push_back(argv[i]);
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
	vector<const eventMetadata*> inputData(inputFiles.size(), 0);
	outputFile->cd();
	for(unsigned int i = 0; i < inputFiles.size(); ++i) {
		inputData[i] = eventMetadata::readEventFile(inputFiles[i]);
		if(not inputData[i]) {
			printErr << "could not read input data from file '" << inputFileNames[i] << "'. Aborting..." << endl;
			return 1;
		}
	}
	eventMetadata* metadata = eventMetadata::merge(inputData, mergeDiffMeta);
	if(not metadata) {
		printErr << "merge failed. Aborting..." << endl;
		return 1;
	}
	outputFile->cd();
	metadata->Write(eventMetadata::objectNameInFile.c_str());
	outputFile->Close();
	for(unsigned int i = 0; i < inputFiles.size(); ++i) {
		inputFiles[i]->Close();
	}
	return 0;
}
