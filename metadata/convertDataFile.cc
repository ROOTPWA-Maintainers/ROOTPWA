
#include <TClonesArray.h>
#include <TFile.h>
#include <TObject.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>

#include "reportingUtilsEnvironment.h"
#include "reportingUtils.hpp"
#include "dataFileWriter.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "convert a data file without meta data to one with meta data" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -i inputFile -o outputFile -l string [-v]" << endl
	     << "    where:" << endl
	     << "        -i file    input file in ROOTPWA format without meta data (default: ./input.root)" << endl
	     << "        -o file    input file in ROOTPWA format with meta data (default: ./output.root)" << endl
	     << "        -l string  label which is saved to the metadata (default: output file name)" << endl
	     << "        -v         verbose; print debug output (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}

int main(int argc, char** argv)
{

// TODO: Add possibility to give a binning map

	printCompilerInfo();
	printLibraryInfo();
	printGitHash();
	cout << endl;

	string prodKinPartNamesObjName = "prodKinParticles";
	string prodKinMomentaLeafName = "prodKinMomenta";
	string decayKinPartNamesObjName = "decayKinParticles";
	string decayKinMomentaLeafName = "decayKinMomenta";
	string inTreeName = "rootPwaEvtTree";

	const string progName = argv[0];
	string inputFileName  = "input.root";
	string outputFileName = "output.root";
	string userString = outputFileName;
	bool debug = false;

	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");

	extern char* optarg;
	int c;
	while((c = getopt(argc, argv, "i:o:l:vh")) != -1)
	{
		switch(c) {
		case 'i':
			inputFileName = optarg;
			break;
		case 'o':
			outputFileName = optarg;
			break;
		case 'l':
			userString = optarg;
			break;
		case 'v':
			debug = true;
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
	TFile* outputFile = TFile::Open(outputFileName.c_str(), "NEW");
	if(not outputFile) {
		printErr << "could not open output file. Aborting..." << endl;
		return 1;
	}
	vector<string> initialStateParticleNames;
	vector<string> finalStateParticleNames;
	{
		TClonesArray* initialStateParticleNames_ = (TClonesArray*)inputFile->Get(prodKinPartNamesObjName.c_str());
		if(not initialStateParticleNames_) {
			printErr << "could not read initial state particle names. Aborting..." << endl;
			return 1;
		}
		TClonesArray* finalStateParticleNames_ = (TClonesArray*)inputFile->Get(decayKinPartNamesObjName.c_str());
		if(not finalStateParticleNames_) {
			printErr << "could not read final state particle names. Aborting..." << endl;
			return 1;
		}
		for(int i = 0; i < initialStateParticleNames_->GetEntries(); ++i) {
			initialStateParticleNames.push_back(string(((TObjString*)(*initialStateParticleNames_)[i])->GetString()));
		}
		if(debug) {
			printDebug << "got initial state particle names " << initialStateParticleNames << endl;
		}
		for(int i = 0; i < finalStateParticleNames_->GetEntries(); ++i) {
			finalStateParticleNames.push_back(string(((TObjString*)(*finalStateParticleNames_)[i])->GetString()));
		}
		if(debug) {
			printDebug << "got final state particle names " << finalStateParticleNames << endl;
		}
	}
	const unsigned int nmbInitialStateParticles = initialStateParticleNames.size();
	const unsigned int nmbFinalStateParticles = finalStateParticleNames.size();
	dataFileWriter fileWriter;
	{
		bool success = fileWriter.initialize(*outputFile,
											 userString,
											 initialStateParticleNames,
											 finalStateParticleNames,
											 map<string, pair<double, double> >(),
											 std::vector<string>(),
											 inTreeName,
											 prodKinMomentaLeafName,
	                                     decayKinMomentaLeafName);
		if(not success) {
			printErr << "could not initialize rootpwaDataFileWriter. Aborting..." << endl;
			return 1;
		}
	}
	TClonesArray* initialStateParticles = 0;
	TClonesArray* finalStateParticles = 0;
	TTree* inputTree = (TTree*)inputFile->Get(inTreeName.c_str());
	inputTree->SetBranchAddress(prodKinMomentaLeafName.c_str(), &initialStateParticles);
	inputTree->SetBranchAddress(decayKinMomentaLeafName.c_str(), &finalStateParticles);
	for(long eventNumber = 0; eventNumber < inputTree->GetEntries(); ++eventNumber) {
		inputTree->GetEntry(eventNumber);
		if(initialStateParticles->GetEntries() != (int)nmbInitialStateParticles) {
			printErr << "received unexpected number of initial state particles (got "
			         << initialStateParticles->GetEntries() << ", expected "
			         << nmbInitialStateParticles << "). Aborting..." << endl;
			return 1;
		}
		if(finalStateParticles->GetEntries() != (int)nmbFinalStateParticles) {
			printErr << "received unexpected number of final state particles (got "
			         << finalStateParticles->GetEntries() << ", expected "
			         << nmbFinalStateParticles << "). Aborting..." << endl;
			return 1;
		}
		vector<TVector3> initialStateMomenta(nmbInitialStateParticles);
		vector<TVector3> finalStateMomenta(nmbFinalStateParticles);
		for(int i = 0; i < initialStateParticles->GetEntries(); ++i) {
			initialStateMomenta[i] = *(TVector3*)(*initialStateParticles)[i];
		}
		for(int i = 0; i < finalStateParticles->GetEntries(); ++i) {
			finalStateMomenta[i] = *(TVector3*)(*finalStateParticles)[i];
		}
		fileWriter.addEvent(initialStateMomenta, finalStateMomenta);
	}
	if(not fileWriter.finalize()) {
		printErr << "finalizing the output file failed." << endl;
		return 1;
	}
	return 0;
}
