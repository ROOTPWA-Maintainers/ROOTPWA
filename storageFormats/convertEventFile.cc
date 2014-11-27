
#include <algorithm>

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TObject.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>

#include "reportingUtilsEnvironment.h"
#include "reportingUtils.hpp"
#include "eventFileWriter.h"


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
	     << "        -t type    type of data (can be 'real', 'generated' or 'accepted', default: 'other')" << endl
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

	string eventsTypeString = "other";

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

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// if the following line is missing, there are error messages of the sort
	// "Warning in <TClass::TClass>: no dictionary for class TVector3 is available",
	// no idea why. Probably there is a correct way of ensuring that the dictionaries
	// are loaded, but I do not know it (yet). KB
	gROOT->ProcessLine("#include <complex>");
#endif

	extern char* optarg;
	int c;
	while((c = getopt(argc, argv, "i:o:l:t:vh")) != -1)
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
		case 't':
			eventsTypeString = optarg;
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
	transform(eventsTypeString.begin(), eventsTypeString.end(), eventsTypeString.begin(), ::tolower);
	eventMetadata::eventsTypeEnum eventsType = eventMetadata::OTHER;
	if(eventsTypeString == "real") {
		eventsType = eventMetadata::REAL;
	} else if(eventsTypeString == "generated") {
		eventsType = eventMetadata::GENERATED;
	} else if(eventsTypeString == "accepted") {
		eventsType = eventMetadata::ACCEPTED;
	} else if(eventsTypeString == "other") {
		// do nothing
	} else {
		printErr << "type '" << eventsTypeString << "' is invalid as an event data type." << endl;
		return 1;
	}

	vector<string> productionKinematicsParticleNames;
	vector<string> decayKinematicsParticleNames;
	{
		TClonesArray* productionKinematicsParticleNames_ = (TClonesArray*)inputFile->Get(prodKinPartNamesObjName.c_str());
		if(not productionKinematicsParticleNames_) {
			printErr << "could not read initial state particle names. Aborting..." << endl;
			return 1;
		}
		TClonesArray* decayKinematicsParticleNames_ = (TClonesArray*)inputFile->Get(decayKinPartNamesObjName.c_str());
		if(not decayKinematicsParticleNames_) {
			printErr << "could not read final state particle names. Aborting..." << endl;
			return 1;
		}
		for(int i = 0; i < productionKinematicsParticleNames_->GetEntries(); ++i) {
			productionKinematicsParticleNames.push_back(string(((TObjString*)(*productionKinematicsParticleNames_)[i])->GetString()));
		}
		if(debug) {
			printDebug << "got initial state particle names " << productionKinematicsParticleNames << endl;
		}
		for(int i = 0; i < decayKinematicsParticleNames_->GetEntries(); ++i) {
			decayKinematicsParticleNames.push_back(string(((TObjString*)(*decayKinematicsParticleNames_)[i])->GetString()));
		}
		if(debug) {
			printDebug << "got final state particle names " << decayKinematicsParticleNames << endl;
		}
	}
	const unsigned int nmbProductionKinematicsParticles = productionKinematicsParticleNames.size();
	const unsigned int nmbDecayKinematicsParticles = decayKinematicsParticleNames.size();
	eventFileWriter fileWriter;
	{
		bool success = fileWriter.initialize(*outputFile,
											 userString,
											 eventsType,
											 productionKinematicsParticleNames,
											 decayKinematicsParticleNames,
											 map<string, pair<double, double> >(),
											 std::vector<string>());
		if(not success) {
			printErr << "could not initialize rootpwaDataFileWriter. Aborting..." << endl;
			return 1;
		}
	}
	TClonesArray* productionKinematicsParticles = 0;
	TClonesArray* decayKinematicsParticles = 0;
	TTree* inputTree = (TTree*)inputFile->Get(inTreeName.c_str());
	inputTree->SetBranchAddress(prodKinMomentaLeafName.c_str(), &productionKinematicsParticles);
	inputTree->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinematicsParticles);
	boost::progress_display progressIndicator(inputTree->GetEntries(), cout, "");
	for(long eventNumber = 0; eventNumber < inputTree->GetEntries(); ++eventNumber) {
		++progressIndicator;
		inputTree->GetEntry(eventNumber);
		if(productionKinematicsParticles->GetEntries() != (int)nmbProductionKinematicsParticles) {
			cout << endl;
			printErr << "received unexpected number of initial state particles (got "
			         << productionKinematicsParticles->GetEntries() << ", expected "
			         << nmbProductionKinematicsParticles << "). Aborting..." << endl;
			return 1;
		}
		if(decayKinematicsParticles->GetEntries() != (int)nmbDecayKinematicsParticles) {
			cout << endl;
			printErr << "received unexpected number of final state particles (got "
			         << decayKinematicsParticles->GetEntries() << ", expected "
			         << nmbDecayKinematicsParticles << "). Aborting..." << endl;
			return 1;
		}
		vector<TVector3> productionKinematicsMomenta(nmbProductionKinematicsParticles);
		vector<TVector3> decayKinematicsMomenta(nmbDecayKinematicsParticles);
		for(int i = 0; i < productionKinematicsParticles->GetEntries(); ++i) {
			productionKinematicsMomenta[i] = *(TVector3*)(*productionKinematicsParticles)[i];
		}
		for(int i = 0; i < decayKinematicsParticles->GetEntries(); ++i) {
			decayKinematicsMomenta[i] = *(TVector3*)(*decayKinematicsParticles)[i];
		}
		fileWriter.addEvent(productionKinematicsMomenta, decayKinematicsMomenta);
	}
	if(not fileWriter.finalize()) {
		printErr << "finalizing the output file failed." << endl;
		return 1;
	}
	inputFile->Close();
	return 0;
}
