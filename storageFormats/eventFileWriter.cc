
#include <map>

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TMD5.h>
#include <TObject.h>
#include <TTree.h>
#include <TVector3.h>

#include "eventFileWriter.h"
#include "eventMetadata.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


void rpwa::md5Wrapper::Update(const double& value) {
	TMD5::Update((UChar_t*)&value, 8);
}

void rpwa::md5Wrapper::Update(const TVector3& vector) {
	Update(vector.X());
	Update(vector.Y());
	Update(vector.Z());
}


rpwa::eventFileWriter::eventFileWriter()
	: _initialized(false),
	  _outfile(0),
	  _eventTree(0),
	  _initialStateMomenta(0),
	  _finalStateMomenta(0),
	  _additionalVariablesToSave(),
	  _metadata(),
	  _nmbInitialStateParticles(0),
	  _nmbFinalStateParticles(0),
	  _md5Calculator() { }


rpwa::eventFileWriter::~eventFileWriter()
{
	reset();
}


bool rpwa::eventFileWriter::initialize(TFile&                                     outputFile,
                                             const string&                              userString,
                                             const vector<string>&                      initialStateParticleNames,
                                             const vector<string>&                      finalStateParticleNames,
                                             const map<string, pair<double, double> >&  binningMap,
                                             const vector<string>&                      additionalVariableLabels,
                                             const string&                              eventTreeName,
                                             const string&                              initialStateMomentaBranchName,
                                             const string&                              finalStateMomentaBranchName,
                                             const string&                              metadataName,
                                             const int&                                 splitlevel,
                                             const int&                                 buffsize)
{
	if(_initialized) {
		printWarn << "trying to initialize when already initialized" << endl;
		return false;
	}
	_outfile = &outputFile;
	_outfile->cd();

	// prepare metadata
	_metadata.setUserString(userString);
	_metadata.setInitialStateParticleNames(initialStateParticleNames);
	_nmbInitialStateParticles = initialStateParticleNames.size();
	_metadata.setFinalStateParticleNames(finalStateParticleNames);
	_nmbFinalStateParticles = finalStateParticleNames.size();
	_metadata.setBinningMap(binningMap);

	// prepare event tree
	_initialStateMomenta = new TClonesArray("TVector3", _nmbInitialStateParticles);
	_finalStateMomenta   = new TClonesArray("TVector3", _nmbFinalStateParticles);
	_eventTree = new TTree(eventTreeName.c_str(), eventTreeName.c_str());
	_eventTree->Branch(initialStateMomentaBranchName.c_str(), "TClonesArray", &_initialStateMomenta, buffsize, splitlevel);
	_eventTree->Branch(finalStateMomentaBranchName.c_str(),   "TClonesArray", &_finalStateMomenta,   buffsize, splitlevel);
	_metadata.setAdditionalSavedVariableLables(additionalVariableLabels);
	_additionalVariablesToSave = vector<double>(additionalVariableLabels.size(), 0.);
	for(unsigned int i = 0; i < additionalVariableLabels.size(); ++i) {
		stringstream strStr;
		strStr << additionalVariableLabels[i] << "/D";
		_eventTree->Branch(additionalVariableLabels[i].c_str(), &_additionalVariablesToSave[i], strStr.str().c_str());
	}

	_initialized = true;
	return _initialized;
}


void rpwa::eventFileWriter::addEvent(const vector<TVector3>& initialStateMomenta,
                                           const vector<TVector3>& finalStateMomenta,
                                           const vector<double>&   additionalVariablesToSave)
{
	if(initialStateMomenta.size() != _nmbInitialStateParticles) {
		printErr << "received unexpected number of initial state particles (got "
		         << initialStateMomenta.size() << ", expected "
		         << _nmbInitialStateParticles << "). Aborting..." << endl;
		throw;
	}
	if(finalStateMomenta.size() != _nmbFinalStateParticles) {
		printErr << "received unexpected number of final state particles (got "
		         << finalStateMomenta.size() << ", expected "
		         << _nmbFinalStateParticles << "). Aborting..." << endl;
		throw;
	}
	if(additionalVariablesToSave.size() != _additionalVariablesToSave.size()) {
		printErr << "received unexpected number of additional variables (got "
		         << additionalVariablesToSave.size() << ", expected "
		         << _additionalVariablesToSave.size() << "). Aborting..." << endl;
		throw;
	}
	for(unsigned int i = 0; i < initialStateMomenta.size(); ++i) {
		const TVector3& initialStateMomentum = initialStateMomenta[i];
		_md5Calculator.Update(initialStateMomentum);
		new ((*_initialStateMomenta)[i]) TVector3(initialStateMomentum);
	}
	for(unsigned int i = 0; i < finalStateMomenta.size(); ++i) {
		const TVector3& finalStateMomentum = finalStateMomenta[i];
		_md5Calculator.Update(finalStateMomentum);
		new ((*_finalStateMomenta)[i]) TVector3(finalStateMomentum);
	}
	for(unsigned int i = 0; i < additionalVariablesToSave.size(); ++i) {
		_md5Calculator.Update(additionalVariablesToSave[i]);
		_additionalVariablesToSave[i] = additionalVariablesToSave[i];
	}
	_eventTree->Fill();
}


bool rpwa::eventFileWriter::finalize() {
	if(not _initialized) {
		printWarn << "trying to finalize when not initialized" << endl;
		return false;
	}
	_outfile->cd();
	_metadata.setContentHash(_md5Calculator.hash());
	_metadata.Write("dataMetadata");
	_outfile->Write();
	_outfile->Close();
	_eventTree = 0;
	reset();
	return true;
}


void rpwa::eventFileWriter::reset() {
	if(_initialStateMomenta) {
		delete _initialStateMomenta;
		_initialStateMomenta = 0;
	}
	if(_finalStateMomenta) {
		delete _finalStateMomenta;
		_finalStateMomenta = 0;
	}
	_eventTree = 0;
	_outfile = 0;
	_initialized = false;
	_metadata = eventMetadata();
}


std::string rpwa::eventFileWriter::calculateHash(TTree* eventTree,
                                                const vector<string> additionalVariableLabels,
                                                const bool&          printProgress,
                                                const string&        initialStateMomentaBranchName,
                                                const string&        finalStateMomentaBranchName)
{
	TClonesArray* initialStateMomenta = 0;
	TClonesArray* finalStateMomenta = 0;
	eventTree->SetBranchAddress(initialStateMomentaBranchName.c_str(), &initialStateMomenta);
	eventTree->SetBranchAddress(finalStateMomentaBranchName.c_str(),   &finalStateMomenta);
	vector<double> additionalVariables(additionalVariableLabels.size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		eventTree->SetBranchAddress(additionalVariableLabels[i].c_str(), &additionalVariables[i]);
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(eventTree->GetEntries(), cout, "") : 0;
	md5Wrapper md5Calculator;
	for(long eventNumber = 0; eventNumber < eventTree->GetEntries(); ++eventNumber) {
		eventTree->GetEntry(eventNumber);
		if(progressIndicator) {
			++(*progressIndicator);
		}
		for(int i = 0; i < initialStateMomenta->GetEntries(); ++i) {
			md5Calculator.Update(*((TVector3*)(*initialStateMomenta)[i]));
		}
		for(int i = 0; i < finalStateMomenta->GetEntries(); ++i) {
			md5Calculator.Update(*((TVector3*)(*finalStateMomenta)[i]));
		}
		for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
			md5Calculator.Update(additionalVariables[i]);
		}
	}
	return md5Calculator.hash();
}
