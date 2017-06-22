
#include "eventMetadata.h"

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include "hashCalculator.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


const std::string rpwa::eventMetadata::objectNameInFile = "eventMetadata";
const std::string rpwa::eventMetadata::eventTreeName = "rootPwaEvtTree";
const std::string rpwa::eventMetadata::productionKinematicsMomentaBranchName = "prodKinMomenta";
const std::string rpwa::eventMetadata::decayKinematicsMomentaBranchName = "decayKinMomenta";


rpwa::eventMetadata::eventMetadata()
	: _userString(""),
	  _contentHash(""),
	  _eventsType(eventMetadata::OTHER),
	  _productionKinematicsParticleNames(),
	  _decayKinematicsParticleNames(),
	  _multibinBoundaries(),
	  _eventTree(0)
{ }


rpwa::eventMetadata::~eventMetadata() { };


ostream& rpwa::eventMetadata::print(ostream& out) const
{
	out << "eventMetadata: " << endl
	    << "    userString ...................... '" << _userString << "'"                  << endl
	    << "    contentHash ..................... '" << _contentHash << "'"                 << endl
	    << "    eventsType ...................... '" << getStringForEventsType(_eventsType) << "'" << endl
	    << "    initial state particle names: ... "  << _productionKinematicsParticleNames  << endl
	    << "    final state particle names: ..... "  << _decayKinematicsParticleNames       << endl
	    << "    multi-bin";
	if(_multibinBoundaries.empty()) {
		out << " ..................... " << "<empty>" << endl;
	} else {
		out << ": " << endl;
		for(multibinBoundariesType::const_iterator it = _multibinBoundaries.begin(); it != _multibinBoundaries.end(); ++it) {
			out << "        variable '" << it->first << "' range " << it->second << endl;
		}
	}
	if(_eventTree) {
		out << "    number of events in file ........ " << _eventTree->GetEntries() << endl;
	}
	out << "    additional branches ............. " << additionalSavedVariableLables() << endl;
	return out;
}


bool rpwa::eventMetadata::operator==(const eventMetadata& rhs) const
{
	if (_userString != rhs._userString) {
		return false;
	}
	if (_contentHash != rhs._contentHash) {
		return false;
	}
	if (_eventsType != rhs._eventsType) {
		return false;
	}
	if (_productionKinematicsParticleNames != rhs._productionKinematicsParticleNames) {
		return false;
	}
	if (_decayKinematicsParticleNames != rhs._decayKinematicsParticleNames) {
		return false;
	}
	if (_multibinBoundaries != rhs._multibinBoundaries) {
		return false;
	}
	if (_additionalSavedVariableLabels != rhs._additionalSavedVariableLabels) {
		return false;
	}
	return true;
}


void rpwa::eventMetadata::appendToUserString(const string& userString,
                                             const string& delimiter)
{
	if(_userString == "") {
		_userString = userString;
	} else {
		stringstream strStr;
		strStr << _userString << delimiter << userString;
		_userString = strStr.str();
	}
}


void rpwa::eventMetadata::setProductionKinematicsParticleNames(const vector<string>& productionKinematicsNames)
{
	_productionKinematicsParticleNames = productionKinematicsNames;
}


void rpwa::eventMetadata::setDecayKinematicsParticleNames(const vector<string>& decayKinematicsParticleNames)
{
	_decayKinematicsParticleNames = decayKinematicsParticleNames;
}


void rpwa::eventMetadata::setBinningVariableLabels(const vector<string>& labels)
{
	for(unsigned int i = 0; i < labels.size(); ++i) {
		_multibinBoundaries[labels[i]] = boundaryType(0., 0.);
	}
}


void rpwa::eventMetadata::setBinningVariableRange(const string& label, const boundaryType& range)
{
	_multibinBoundaries[label] = range;
}


void rpwa::eventMetadata::setMultibinBoundaries(const multibinBoundariesType& multibinBoundaries)
{
	_multibinBoundaries = multibinBoundaries;
}


string rpwa::eventMetadata::recalculateHash(const bool& printProgress) const
{
	TClonesArray* productionKinematicsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	hashCalculator hashor;
	if(not _eventTree) {
		printWarn << "input tree not found in metadata." << endl;
		return "";
	}
	if(_eventTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0)
	{
		printWarn << "could not set address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
		return "";
	}
	if(_eventTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta)) {
		printWarn << "could not set address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
		return "";
	}
	vector<double> additionalVariables(additionalSavedVariableLables().size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		if(_eventTree->SetBranchAddress(additionalSavedVariableLables()[i].c_str(), &additionalVariables[i]) < 0) {
			printWarn << "could not set address for branch '" << additionalSavedVariableLables()[i].c_str() << "'." << endl;
			return "";
		}
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(_eventTree->GetEntries(), cout, "") : 0;
	for(long eventNumber = 0; eventNumber < _eventTree->GetEntries(); ++eventNumber) {
		_eventTree->GetEntry(eventNumber);
		if(progressIndicator) {
			++(*progressIndicator);
		}
		for(int i = 0; i < productionKinematicsMomenta->GetEntries(); ++i) {
			hashor.Update(*((TVector3*)(*productionKinematicsMomenta)[i]));
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


Long64_t rpwa::eventMetadata::Merge(TCollection* /*list*/, Option_t* /*option*/) {
	printErr << "data files cannot be merged with hadd. Please use $ROOTPWA/build/bin/mergeDatafiles." << endl;
	throw;
}


eventMetadata* rpwa::eventMetadata::merge(const vector<const eventMetadata*>& inputData,
                                          const bool mergeDiffMeta,
                                          const int& splitlevel,
                                          const int& buffsize)
{
	eventMetadata* mergee = new eventMetadata();
	if(inputData.empty()) {
		printWarn << "trying to merge without input data." << endl;
		return 0;
	}
	hashCalculator hashor;
	const unsigned int nmbProductionKinematicsParticles = inputData[0]->productionKinematicsParticleNames().size();
	const unsigned int nmbDecayKinematicsParticles = inputData[0]->decayKinematicsParticleNames().size();
	TClonesArray* productionKinematicsMomenta = new TClonesArray("TVector3", nmbProductionKinematicsParticles);
	TClonesArray* decayKinematicsMomenta   = new TClonesArray("TVector3", nmbDecayKinematicsParticles);
	mergee->_eventTree = new TTree(eventTreeName.c_str(), eventTreeName.c_str());
	mergee->_eventTree->Branch(eventMetadata::productionKinematicsMomentaBranchName.c_str(), "TClonesArray", &productionKinematicsMomenta, buffsize, splitlevel);
	mergee->_eventTree->Branch(eventMetadata::decayKinematicsMomentaBranchName.c_str(),   "TClonesArray", &decayKinematicsMomenta,   buffsize, splitlevel);
	vector<double> additionalSavedVariables;
	bool first = true;
	multibinBoundariesType mergedmultibinBoundaries;
	for(unsigned int inputDataNumber = 0; inputDataNumber < inputData.size(); ++inputDataNumber) {
		const eventMetadata* metadata = inputData[inputDataNumber];
		TTree* inputTree = metadata->eventTree();
		if(not inputTree) {
			printWarn << "got NULL-pointer to inputTree when merging." << endl;
			delete mergee->_eventTree;
			delete mergee;
			return 0;
		}
		if(first) {
			first = false;
			if((mergee->productionKinematicsParticleNames().empty()) and
			   (mergee->decayKinematicsParticleNames().empty()))
			{
				mergee->setProductionKinematicsParticleNames(metadata->productionKinematicsParticleNames());
				mergee->setDecayKinematicsParticleNames(metadata->decayKinematicsParticleNames());
			}
			if(not mergeDiffMeta) {
				if(mergee->multibinBoundaries().empty()) {
					mergee->setMultibinBoundaries(metadata->multibinBoundaries());
				}
			} else {
				mergedmultibinBoundaries = metadata->multibinBoundaries();
			}
			if(mergee->additionalSavedVariableLables().empty()) {
				mergee->setAdditionalSavedVariableLables(metadata->additionalSavedVariableLables());
				additionalSavedVariables.resize(mergee->additionalSavedVariableLables().size(), 0.);
				for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
					stringstream strStr;
					strStr << mergee->additionalSavedVariableLables()[i] << "/D";
					mergee->_eventTree->Branch(mergee->additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i], strStr.str().c_str());
				}
			}
		}
		mergee->appendToUserString(metadata->userString());
		if(mergee->productionKinematicsParticleNames() != metadata->productionKinematicsParticleNames()) {
			printWarn << "particle names of production kinematics differ." << endl;
			delete mergee->_eventTree;
			delete mergee;
			return 0;
		}
		if(mergee->decayKinematicsParticleNames() != metadata->decayKinematicsParticleNames()) {
			printWarn << "particle names of decay kinematics differ." << endl;
			delete mergee->_eventTree;
			delete mergee;
			return 0;
		}
		if (not mergeDiffMeta) {
			if(mergee->multibinBoundaries() != metadata->multibinBoundaries()) {
				printWarn << "multibin ranges differ." << endl;
				delete mergee->_eventTree;
				delete mergee;
				return 0;
			}
		} else {
			for(multibinBoundariesType::const_iterator iterator = mergedmultibinBoundaries.begin(); iterator != mergedmultibinBoundaries.end(); iterator++) {
				const string& binningVariable = iterator->first;
				const boundaryType& binningRange = iterator->second;
				const multibinBoundariesType::const_iterator toAdd = metadata->multibinBoundaries().find(binningVariable);
				if (toAdd != metadata->multibinBoundaries().end()) { // check if binningVariable exists (uniquely) in other multibinBoundaries
					if (binningRange.first > toAdd->second.first) {
						mergedmultibinBoundaries[binningVariable].first = toAdd->second.first;
					}
					if (binningRange.second < toAdd->second.second) {
						mergedmultibinBoundaries[binningVariable].second = toAdd->second.second;
					}
				} else {
					printWarn << "files do not use the same binning variables." << endl;
					delete mergee->_eventTree;
					delete mergee;
					return 0;
				}
			}
		}
		if(inputTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
			delete mergee->_eventTree;
			delete mergee;
			return 0;
		}
		if(inputTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
			delete mergee->_eventTree;
			delete mergee;
			return 0;
		}
		for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
			if(inputTree->SetBranchAddress(mergee->additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i]) < 0)
			{
				printWarn << "could not set address for branch '" << mergee->additionalSavedVariableLables()[i] << "'." << endl;
				delete mergee->_eventTree;
				delete mergee;
				return 0;
			}
		}
		for(long eventNumber = 0; eventNumber < inputTree->GetEntries(); ++eventNumber) {
			inputTree->GetEntry(eventNumber);
			for(int i = 0; i < productionKinematicsMomenta->GetEntries(); ++i) {
				hashor.Update(*((TVector3*)(*productionKinematicsMomenta)[i]));
			}
			for(int i = 0; i < decayKinematicsMomenta->GetEntries(); ++i) {
				hashor.Update(*((TVector3*)(*decayKinematicsMomenta)[i]));
			}
			for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
				hashor.Update(additionalSavedVariables[i]);
			}
			mergee->_eventTree->Fill();
		}
	}
	if(mergeDiffMeta) {
		mergee->setMultibinBoundaries(mergedmultibinBoundaries);
	}
	mergee->setContentHash(hashor.hash());
	return mergee;
}


const eventMetadata* rpwa::eventMetadata::readEventFile(TFile* inputFile, const bool& quiet)
{
	eventMetadata* eventMeta = (eventMetadata*)inputFile->Get(objectNameInFile.c_str());
	if(not eventMeta) {
		if(not quiet) {
			printWarn << "could not find event metadata." << endl;
		}
		return 0;
	}
	eventMeta->_eventTree = (TTree*)inputFile->Get(eventTreeName.c_str());
	if(not eventMeta->_eventTree) {
		if(not quiet) {
			printWarn << "could not find event tree." << endl;
		}
		return 0;
	}
	return eventMeta;
}


Int_t rpwa::eventMetadata::Write(const char* name, Int_t option, Int_t bufsize) const
{
	Int_t retval = 0;
	if(_eventTree) {
		retval = _eventTree->Write();
	}
	return retval + TObject::Write(name, option, bufsize);
}


std::string rpwa::eventMetadata::getStringForEventsType(const eventsTypeEnum& type)
{
	switch(type) {
		case OTHER:
			return "other";
		case REAL:
			return "real";
		case GENERATED:
			return "generated";
		case ACCEPTED:
			return "accepted";
	}
	return "UNKNOWN";
}
