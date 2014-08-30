
#include <map>

#include <TClonesArray.h>
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <TVector3.h>

#include "dataFileWriter.h"
#include "dataMetadata.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


rpwa::rootpwaDataFileWriter::rootpwaDataFileWriter()
	: _initialized(false),
	  _outfile(0),
	  _eventTree(0),
	  _initialStateMomenta(0),
	  _finalStateMomenta(0),
	  _additionalVariablesToSave(),
	  _metadata(),
	  _nmbInitialStateParticles(0),
	  _nmbFinalStateParticles(0) { }


rpwa::rootpwaDataFileWriter::~rootpwaDataFileWriter()
{
	reset();
}


bool rpwa::rootpwaDataFileWriter::initialize(TFile&                                     outputFile,
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
	_additionalVariablesToSave = vector<double>(additionalVariableLabels.size(), 0.);
	for(unsigned int i = 0; i < additionalVariableLabels.size(); ++i) {
		stringstream strStr;
		strStr << additionalVariableLabels[i] << "/D";
		_eventTree->Branch(additionalVariableLabels[i].c_str(), &_additionalVariablesToSave[i], strStr.str().c_str());
	}

	_initialized = true;
	return _initialized;
}


void rpwa::rootpwaDataFileWriter::addEvent(const vector<TVector3>& initialStateMomenta,
                                           const vector<TVector3>& finalStateMomenta,
                                           const vector<double>&   additionalVariablesToSave)
{
	// TODO: calculate hash
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
		new ((*_initialStateMomenta)[i]) TVector3(initialStateMomenta[i]);
	}
	for(unsigned int i = 0; i < finalStateMomenta.size(); ++i) {
		new ((*_finalStateMomenta)[i]) TVector3(finalStateMomenta[i]);
	}
	for(unsigned int i = 0; i < additionalVariablesToSave.size(); ++i) {
		_additionalVariablesToSave[i] = additionalVariablesToSave[i];
	}
	_eventTree->Fill();
}


bool rpwa::rootpwaDataFileWriter::finalize() {
	if(not _initialized) {
		printWarn << "trying to finalize when not initialized" << endl;
		return false;
	}
	_outfile->cd();
	_metadata.Write("dataMetadata");
	_outfile->Write();
	_outfile->Close();
	_eventTree = 0;
	reset();
	return true;
}


void rpwa::rootpwaDataFileWriter::reset() {
	if(_eventTree) {
		delete _eventTree;
		_eventTree = 0;
	}
	if(_initialStateMomenta) {
		delete _initialStateMomenta;
		_initialStateMomenta = 0;
	}
	if(_finalStateMomenta) {
		delete _finalStateMomenta;
		_finalStateMomenta = 0;
	}
	_outfile = 0;
	_initialized = false;
	_metadata = dataMetadata();
}
