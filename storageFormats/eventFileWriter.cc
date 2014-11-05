
#include <map>

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
	  _outputFile(0),
	  _metadata(),
	  _productionKinematicsMomenta(0),
	  _decayKinematicsMomenta(0),
	  _additionalVariablesToSave(),
	  _nmbProductionKinematicsParticles(0),
	  _nmbDecayKinematicsParticles(0),
	  _hashCalculator() { }


rpwa::eventFileWriter::~eventFileWriter()
{
	reset();
}


bool rpwa::eventFileWriter::initialize(TFile&                                     outputFile,
                                       const string&                              userString,
                                       const eventMetadata::eventsTypeEnum&       eventsType,
                                       const vector<string>&                      productionKinematicsParticleNames,
                                       const vector<string>&                      decayKinematicsParticleNames,
                                       const map<string, pair<double, double> >&  binningMap,
                                       const vector<string>&                      additionalVariableLabels,
                                       const int&                                 splitlevel,
                                       const int&                                 buffsize)
{
	if(_initialized) {
		printWarn << "trying to initialize when already initialized." << endl;
		return false;
	}
	_outputFile = &outputFile;
	_outputFile->cd();

	// prepare metadata
	_metadata.setUserString(userString);
	_metadata.setEventsType(eventsType);
	_metadata.setProductionKinematicsParticleNames(productionKinematicsParticleNames);
	_nmbProductionKinematicsParticles = productionKinematicsParticleNames.size();
	_metadata.setDecayKinematicsParticleNames(decayKinematicsParticleNames);
	_nmbDecayKinematicsParticles = decayKinematicsParticleNames.size();
	_metadata.setBinningMap(binningMap);

	// prepare event tree
	_productionKinematicsMomenta = new TClonesArray("TVector3", _nmbProductionKinematicsParticles);
	_decayKinematicsMomenta   = new TClonesArray("TVector3", _nmbDecayKinematicsParticles);
	_metadata._eventTree = new TTree(eventMetadata::eventTreeName.c_str(), eventMetadata::eventTreeName.c_str());
	_metadata._eventTree->Branch(eventMetadata::productionKinematicsMomentaBranchName.c_str(), "TClonesArray", &_productionKinematicsMomenta, buffsize, splitlevel);
	_metadata._eventTree->Branch(eventMetadata::decayKinematicsMomentaBranchName.c_str(),   "TClonesArray", &_decayKinematicsMomenta,   buffsize, splitlevel);
	_metadata.setAdditionalSavedVariableLables(additionalVariableLabels);
	_additionalVariablesToSave = vector<double>(additionalVariableLabels.size(), 0.);
	for(unsigned int i = 0; i < additionalVariableLabels.size(); ++i) {
		stringstream strStr;
		strStr << additionalVariableLabels[i] << "/D";
		_metadata._eventTree->Branch(additionalVariableLabels[i].c_str(), &_additionalVariablesToSave[i], strStr.str().c_str());
	}

	_initialized = true;
	return _initialized;
}


void rpwa::eventFileWriter::addEvent(const TClonesArray&   productionKinematicsMomenta,
                                     const TClonesArray&   decayKinematicsMomenta,
                                     const vector<double>& additionalVariablesToSave)
{
	if(not _initialized) {
		printWarn << "trying to add event when not initialized." << endl;
		return;
	}
	if(productionKinematicsMomenta.GetEntries() != (int) _nmbProductionKinematicsParticles) {
		printErr << "received unexpected number of initial state particles (got "
		         << productionKinematicsMomenta.GetEntries() << ", expected "
		         << _nmbProductionKinematicsParticles << "). Aborting..." << endl;
		throw;
	}
	if(decayKinematicsMomenta.GetEntries() != (int) _nmbDecayKinematicsParticles) {
		printErr << "received unexpected number of final state particles (got "
		         << decayKinematicsMomenta.GetEntries() << ", expected "
		         << _nmbDecayKinematicsParticles << "). Aborting..." << endl;
		throw;
	}
	if(additionalVariablesToSave.size() != _additionalVariablesToSave.size()) {
		printErr << "received unexpected number of additional variables (got "
		         << additionalVariablesToSave.size() << ", expected "
		         << _additionalVariablesToSave.size() << "). Aborting..." << endl;
		throw;
	}
	for(int i = 0; i < productionKinematicsMomenta.GetEntries(); ++i) {
		const TVector3& productionKinematicsMomentum = *((TVector3*) productionKinematicsMomenta[i]);
		_hashCalculator.Update(productionKinematicsMomentum);
		new ((*_productionKinematicsMomenta)[i]) TVector3(productionKinematicsMomentum);
	}
	for(int i = 0; i < decayKinematicsMomenta.GetEntries(); ++i) {
		const TVector3& decayKinematicsMomentum = *((TVector3*) decayKinematicsMomenta[i]);
		_hashCalculator.Update(decayKinematicsMomentum);
		new ((*_decayKinematicsMomenta)[i]) TVector3(decayKinematicsMomentum);
	}
	for(unsigned int i = 0; i < additionalVariablesToSave.size(); ++i) {
		_hashCalculator.Update(additionalVariablesToSave[i]);
		_additionalVariablesToSave[i] = additionalVariablesToSave[i];
	}
	_metadata._eventTree->Fill();
}


void rpwa::eventFileWriter::addEvent(const vector<TVector3>& productionKinematicsMomenta,
                                     const vector<TVector3>& decayKinematicsMomenta,
                                     const vector<double>&   additionalVariablesToSave)
{
	TClonesArray prodKinClonesArray = TClonesArray("TVector3", productionKinematicsMomenta.size());
	TClonesArray decayKinClonesArray = TClonesArray("TVector3", decayKinematicsMomenta.size());

	for(unsigned int i = 0; i < productionKinematicsMomenta.size(); ++i) {
		new(prodKinClonesArray[i]) TVector3(productionKinematicsMomenta[i]);
	}
	for(unsigned int i = 0; i < decayKinematicsMomenta.size(); ++i) {
		new(decayKinClonesArray[i]) TVector3(decayKinematicsMomenta[i]);
	}
	this->addEvent(prodKinClonesArray, decayKinClonesArray, additionalVariablesToSave);
}


bool rpwa::eventFileWriter::finalize() {
	if(not _initialized) {
		printWarn << "trying to finalize when not initialized." << endl;
		return false;
	}
	_metadata.setContentHash(_hashCalculator.hash());
	_outputFile->cd();
	_metadata.Write(eventMetadata::objectNameInFile.c_str());
	_outputFile->Close();
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
	_outputFile = 0;
	_hashCalculator = hashCalculator();
	_initialized = false;
}
