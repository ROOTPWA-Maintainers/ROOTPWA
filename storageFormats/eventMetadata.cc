
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
	: _auxString(""),
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
	    << "    auxString ....................... '" << _auxString << "'"                   << endl
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
	out << "    additional branches ............. " << additionalTreeVariableNames() << endl;
	out << "    auxValues";
	if(_auxValues.empty()) {
		out << " ....................... " << "<empty>" << endl;
	} else {
		out << ": " << endl;
		for(const auto& auxValue: _auxValues) {
			out << "        variable '" << auxValue.first << "' value " << auxValue.second << endl;
		}
		cout << endl;
	}
	return out;
}


bool rpwa::eventMetadata::operator==(const eventMetadata& rhs) const
{
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
	if (_additionalTreeVariableNames != rhs._additionalTreeVariableNames) {
		return false;
	}
	return true;
}


void rpwa::eventMetadata::appendToAuxString(const string& auxString,
                                            const string& delimiter)
{
	if(_auxString == "") {
		_auxString = auxString;
	} else {
		stringstream strStr;
		strStr << _auxString << delimiter << auxString;
		_auxString = strStr.str();
	}
}


void rpwa::eventMetadata::setBinningVariableLabels(const vector<string>& labels)
{
	for(unsigned int i = 0; i < labels.size(); ++i) {
		_multibinBoundaries[labels[i]] = boundaryType(0., 0.);
	}
}


bool rpwa::eventMetadata::updateHashor(hashCalculator& hashor, const bool& printProgress) const
{
	TClonesArray* productionKinematicsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	if(not _eventTree) {
		printWarn << "input tree not found in metadata." << endl;
		return false;
	}
	if(_eventTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0)
	{
		printWarn << "could not set address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
		return false;
	}
	if(_eventTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta)) {
		printWarn << "could not set address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
		return false;
	}
	vector<double> additionalVariables(additionalTreeVariableNames().size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		if(_eventTree->SetBranchAddress(additionalTreeVariableNames()[i].c_str(), &additionalVariables[i]) < 0) {
			printWarn << "could not set address for branch '" << additionalTreeVariableNames()[i].c_str() << "'." << endl;
			return false;
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
	return true;
}


string rpwa::eventMetadata::recalculateHash(const bool& printProgress) const {
	hashCalculator hashor;
	if (updateHashor(hashor, printProgress)) {
		return hashor.hash();
	} else {
		return "";
	}
}


Long64_t rpwa::eventMetadata::Merge(TCollection* /*list*/, Option_t* /*option*/) {
	printErr << "data files cannot be merged with hadd. Please use $ROOTPWA/build/bin/mergeDatafiles." << endl;
	throw;
}


eventMetadata* rpwa::eventMetadata::merge(const vector<const eventMetadata*>& inputData,
                                          const bool mergeBinBoundaries,
                                          const bool mergeAuxString,
                                          const bool mergeAuxValues,
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
	vector<double> additionalTreeVariables;
	bool first = true;
	multibinBoundariesType mergedMultibinBoundaries;
	for(unsigned int inputDataNumber = 0; inputDataNumber < inputData.size(); ++inputDataNumber) {
		const eventMetadata* metadata = inputData[inputDataNumber];
		TTree* inputTree = metadata->eventTree();
		if(not inputTree) {
			printWarn << "got NULL-pointer to inputTree when merging." << endl;
			goto mergeFailed;
		}
		if(first) {
			first = false;
			mergee->setProductionKinematicsParticleNames(metadata->productionKinematicsParticleNames());
			mergee->setDecayKinematicsParticleNames(metadata->decayKinematicsParticleNames());
			if(not mergeBinBoundaries) {
				mergee->setMultibinBoundaries(metadata->multibinBoundaries());
			} else {
				mergedMultibinBoundaries = metadata->multibinBoundaries();
			}
			mergee->setAdditionalTreeVariableNames(metadata->additionalTreeVariableNames());
			additionalTreeVariables.resize(mergee->additionalTreeVariableNames().size(), 0.);
			for(unsigned int i = 0; i < additionalTreeVariables.size(); ++i) {
				stringstream strStr;
				strStr << mergee->additionalTreeVariableNames()[i] << "/D";
				mergee->_eventTree->Branch(mergee->additionalTreeVariableNames()[i].c_str(), &additionalTreeVariables[i], strStr.str().c_str());
			}
		}
		if (mergeAuxString)
			mergee->appendToAuxString(metadata->auxString());
		if(mergee->productionKinematicsParticleNames() != metadata->productionKinematicsParticleNames()) {
			printWarn << "particle names of production kinematics differ." << endl;
			goto mergeFailed;
		}
		if(mergee->decayKinematicsParticleNames() != metadata->decayKinematicsParticleNames()) {
			printWarn << "particle names of decay kinematics differ." << endl;
			goto mergeFailed;
		}
		if (not mergeBinBoundaries) {
			if(mergee->multibinBoundaries() != metadata->multibinBoundaries()) {
				printWarn << "multibin ranges differ." << endl;
				goto mergeFailed;
			}
		} else {
			for(multibinBoundariesType::const_iterator iterator = mergedMultibinBoundaries.begin(); iterator != mergedMultibinBoundaries.end(); iterator++) {
				const string& binningVariable = iterator->first;
				const boundaryType& binningRange = iterator->second;
				const multibinBoundariesType::const_iterator toAdd = metadata->multibinBoundaries().find(binningVariable);
				if (toAdd != metadata->multibinBoundaries().end()) { // check if binningVariable exists (uniquely) in other multibinBoundaries
					if (toAdd->second.first < binningRange.first) {
						mergedMultibinBoundaries[binningVariable].first = toAdd->second.first;
					}
					if (toAdd->second.second > binningRange.second) {
						mergedMultibinBoundaries[binningVariable].second = toAdd->second.second;
					}
				} else {
					printWarn << "files do not use the same binning variables." << endl;
					goto mergeFailed;
				}
			}
		}
		if(inputTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
			goto mergeFailed;
		}
		if(inputTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
			goto mergeFailed;
		}
		for(unsigned int i = 0; i < additionalTreeVariables.size(); ++i) {
			if(inputTree->SetBranchAddress(mergee->additionalTreeVariableNames()[i].c_str(), &additionalTreeVariables[i]) < 0)
			{
				printWarn << "could not set address for branch '" << mergee->additionalTreeVariableNames()[i] << "'." << endl;
				goto mergeFailed;
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
			for(unsigned int i = 0; i < additionalTreeVariables.size(); ++i) {
				hashor.Update(additionalTreeVariables[i]);
			}
			mergee->_eventTree->Fill();
		}

		if (mergeAuxValues) {
			for (const auto& nameValue : metadata->auxValues()) {
				if (not mergee->hasAuxValue(nameValue.first)) {
					mergee->setAuxValue(nameValue.first, nameValue.second);
				} else { // the auxiliary value exists in multiple files -> only merge if it is the same
					if (mergee->auxValue(nameValue.first) != nameValue.second) {
						printWarn << "auxiliary value '" << nameValue.first << "' has different values in the different files." << endl;
						goto mergeFailed;
					}
				}
			}
		}
	}

	if(mergeBinBoundaries) {
		mergee->setMultibinBoundaries(mergedMultibinBoundaries);
	}

	mergee->setContentHash(hashor.hash());
	return mergee;
mergeFailed:
	delete mergee->_eventTree;
	delete mergee;
	return 0;
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
	return "unknown";
}


bool
additionalTreeVariables::setBranchAddresses(const eventMetadata& metaData)
{
	TTree* tree = metaData.eventTree();

	if (tree != nullptr) {
		for (const auto& name: metaData.additionalTreeVariableNames()) {
			int err = tree->SetBranchAddress(name.c_str(), &(_additionalTreeVariables[name]));
			if (err < 0){
				printErr << "could not set branch address for branch '" << name << "' (error code " << err << ")." << endl;
				_additionalTreeVariables.clear();
				return false;
			}
		}
	} else {
		printErr << "no tree in eventMetadata object." << std::endl;
		return false;
	}

	return true;
}


bool
additionalTreeVariables::inBoundaries(const multibinBoundariesType& boundaries) const
{
	for (const auto& elem: boundaries) {
		const std::string& name = elem.first;
		const boundaryType& range = elem.second;
		const std::map<std::string, double>::const_iterator it = _additionalTreeVariables.find(name);
		if (it != _additionalTreeVariables.end()) {
			const double variable = it->second;
			if (variable < range.first or variable >= range.second) {
				return false;
			}
		} else {
			printErr << "binning variable '" << name << "' not in additional tree variables." << std::endl;
			return false;
		}
	}

	return true;
}
