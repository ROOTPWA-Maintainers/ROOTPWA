
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


rpwa::eventFileWriter::eventFileWriter()
	: _initialized(false),
	  _outfile(0),
	  _eventTree(0),
	  _productionKinematicsMomenta(0),
	  _decayKinematicsMomenta(0),
	  _additionalVariablesToSave(),
	  _metadata(),
	  _nmbProductionKinematicsParticles(0),
	  _nmbDecayKinematicsParticles(0),
	  _hashCalculator() { }


rpwa::eventFileWriter::~eventFileWriter()
{
	reset();
}


bool rpwa::eventFileWriter::initialize(TFile&                                     outputFile,
                                       const string&                              userString,
                                       const vector<string>&                      productionKinematicsParticleNames,
                                       const vector<string>&                      decayKinematicsParticleNames,
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
	_metadata.setProductionKinematicsParticleNames(productionKinematicsParticleNames);
	_nmbProductionKinematicsParticles = productionKinematicsParticleNames.size();
	_metadata.setDecayKinematicsParticleNames(decayKinematicsParticleNames);
	_nmbDecayKinematicsParticles = decayKinematicsParticleNames.size();
	_metadata.setBinningMap(binningMap);

	// prepare event tree
	_productionKinematicsMomenta = new TClonesArray("TVector3", _nmbProductionKinematicsParticles);
	_decayKinematicsMomenta   = new TClonesArray("TVector3", _nmbDecayKinematicsParticles);
	_eventTree = new TTree(eventTreeName.c_str(), eventTreeName.c_str());
	_eventTree->Branch(initialStateMomentaBranchName.c_str(), "TClonesArray", &_productionKinematicsMomenta, buffsize, splitlevel);
	_eventTree->Branch(finalStateMomentaBranchName.c_str(),   "TClonesArray", &_decayKinematicsMomenta,   buffsize, splitlevel);
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


void rpwa::eventFileWriter::addEvent(const vector<TVector3>&       productionKinematicsMomenta,
                                           const vector<TVector3>& decayKinematicsMomenta,
                                           const vector<double>&   additionalVariablesToSave)
{
	if(productionKinematicsMomenta.size() != _nmbProductionKinematicsParticles) {
		printErr << "received unexpected number of initial state particles (got "
		         << productionKinematicsMomenta.size() << ", expected "
		         << _nmbProductionKinematicsParticles << "). Aborting..." << endl;
		throw;
	}
	if(decayKinematicsMomenta.size() != _nmbDecayKinematicsParticles) {
		printErr << "received unexpected number of final state particles (got "
		         << decayKinematicsMomenta.size() << ", expected "
		         << _nmbDecayKinematicsParticles << "). Aborting..." << endl;
		throw;
	}
	if(additionalVariablesToSave.size() != _additionalVariablesToSave.size()) {
		printErr << "received unexpected number of additional variables (got "
		         << additionalVariablesToSave.size() << ", expected "
		         << _additionalVariablesToSave.size() << "). Aborting..." << endl;
		throw;
	}
	for(unsigned int i = 0; i < productionKinematicsMomenta.size(); ++i) {
		const TVector3& productionKinematicsMomentum = productionKinematicsMomenta[i];
		_hashCalculator.Update(productionKinematicsMomentum);
		new ((*_productionKinematicsMomenta)[i]) TVector3(productionKinematicsMomentum);
	}
	for(unsigned int i = 0; i < decayKinematicsMomenta.size(); ++i) {
		const TVector3& decayKinematicsMomentum = decayKinematicsMomenta[i];
		_hashCalculator.Update(decayKinematicsMomentum);
		new ((*_decayKinematicsMomenta)[i]) TVector3(decayKinematicsMomentum);
	}
	for(unsigned int i = 0; i < additionalVariablesToSave.size(); ++i) {
		_hashCalculator.Update(additionalVariablesToSave[i]);
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
	_metadata.setContentHash(_hashCalculator.hash());
	_metadata.Write("dataMetadata");
	_outfile->Write();
	_outfile->Close();
	_eventTree = 0;
	reset();
	return true;
}


void rpwa::eventFileWriter::reset() {
	if(_productionKinematicsMomenta) {
		delete _productionKinematicsMomenta;
		_productionKinematicsMomenta = 0;
	}
	if(_decayKinematicsMomenta) {
		delete _decayKinematicsMomenta;
		_decayKinematicsMomenta = 0;
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
	TClonesArray* productionKinameticsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	eventTree->SetBranchAddress(initialStateMomentaBranchName.c_str(), &productionKinameticsMomenta);
	eventTree->SetBranchAddress(finalStateMomentaBranchName.c_str(),   &decayKinematicsMomenta);
	vector<double> additionalVariables(additionalVariableLabels.size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		eventTree->SetBranchAddress(additionalVariableLabels[i].c_str(), &additionalVariables[i]);
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(eventTree->GetEntries(), cout, "") : 0;
	hashCalculator hashor;
	for(long eventNumber = 0; eventNumber < eventTree->GetEntries(); ++eventNumber) {
		eventTree->GetEntry(eventNumber);
		if(progressIndicator) {
			++(*progressIndicator);
		}
		for(int i = 0; i < productionKinameticsMomenta->GetEntries(); ++i) {
			hashor.Update(*((TVector3*)(*productionKinameticsMomenta)[i]));
		}
		for(int i = 0; i < decayKinematicsMomenta->GetEntries(); ++i) {
			hashor.Update(*((TVector3*)(*decayKinematicsMomenta)[i]));
		}
		for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
			hashor.Update(additionalVariables[i]);
		}
	}
	return hashor.hash();
}
