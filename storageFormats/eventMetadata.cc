
#include "eventMetadata.h"

#include <boost/progress.hpp>

#include <TClonesArray.h>
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
	  _productionKinematicsParticleNames(),
	  _decayKinematicsParticleNames(),
	  _binningMap()
{ }


rpwa::eventMetadata::~eventMetadata() { };


ostream& rpwa::eventMetadata::print(ostream& out) const
{
	out << "dataMetadata: " << endl
	    << "    userString ...................... '" << _userString << "'"         << endl
	    << "    contentHash ..................... '" << _contentHash << "'"        << endl
	    << "    initial state particle names: ... "  << _productionKinematicsParticleNames << endl
        << "    final state particle names: ..... "  << _decayKinematicsParticleNames   << endl
	    << "    binning map";
	if(_binningMap.empty()) {
		out << " ..................... " << "<empty>" << endl;
	} else {
		out << ": " << endl;
		for(binningMapType::const_iterator it = _binningMap.begin(); it != _binningMap.end(); ++it) {
			out << "        variable '" << it->first << "' range " << it->second << endl;
		}
	}
	return out;
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
		_binningMap[labels[i]] = rangePairType(0., 0.);
	}
}


void rpwa::eventMetadata::setBinningVariableRange(const string& label, const rangePairType& range)
{
	_binningMap[label] = range;
}


void rpwa::eventMetadata::setBinningMap(const binningMapType& binningMap)
{
	_binningMap = binningMap;
}


string rpwa::eventMetadata::hash(TTree* eventTree, const bool& printProgress) const
{
	TClonesArray* productionKinematicsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	hashCalculator hashor;
	if(eventTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0)
	{
		printWarn << "could not set address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
		return "";
	}
	if(eventTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta)) {
		printWarn << "could not set address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
		return "";
	}
	vector<double> additionalVariables(additionalSavedVariableLables().size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		if(eventTree->SetBranchAddress(additionalSavedVariableLables()[i].c_str(), &additionalVariables[i]) < 0) {
			printWarn << "could not set address for branch '" << additionalSavedVariableLables()[i].c_str() << "'." << endl;
			return "";
		}
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(eventTree->GetEntries(), cout, "") : 0;
	for(long eventNumber = 0; eventNumber < eventTree->GetEntries(); ++eventNumber) {
		eventTree->GetEntry(eventNumber);
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


Long64_t rpwa::eventMetadata::Merge(TCollection* list, Option_t* option) {
	printErr << "data files cannot be merged with hadd. Please use $ROOTPWA/build/bin/mergeDatafiles." << endl;
	throw;
}


TTree* rpwa::eventMetadata::merge(const vector<pair<const eventMetadata*, TTree*> >& inputData,
                                  const int& splitlevel,
                                  const int& buffsize)
{
	if(inputData.empty()) {
		printWarn << "trying to merge without input data." << endl;
		return false;
	}
	hashCalculator hashor;
	const unsigned int nmbProductionKinematicsParticles = inputData[0].first->productionKinematicsParticleNames().size();
	const unsigned int nmbDecayKinematicsParticles = inputData[0].first->decayKinematicsParticleNames().size();
	TClonesArray* productionKinematicsMomenta = new TClonesArray("TVector3", nmbProductionKinematicsParticles);
	TClonesArray* decayKinematicsMomenta   = new TClonesArray("TVector3", nmbDecayKinematicsParticles);
	TTree* outputTree = new TTree(eventTreeName.c_str(), eventTreeName.c_str());
	outputTree->Branch(eventMetadata::productionKinematicsMomentaBranchName.c_str(), "TClonesArray", &productionKinematicsMomenta, buffsize, splitlevel);
	outputTree->Branch(eventMetadata::decayKinematicsMomentaBranchName.c_str(),   "TClonesArray", &decayKinematicsMomenta,   buffsize, splitlevel);
	vector<double> additionalSavedVariables;
	bool first = true;
	for(unsigned int inputDataNumber = 0; inputDataNumber < inputData.size(); ++inputDataNumber) {
		const eventMetadata* metadata = inputData[inputDataNumber].first;
		TTree* inputTree = inputData[inputDataNumber].second;
		if(not inputTree) {
			printWarn << "got NULL-pointer to inputTree when merging." << endl;
			delete outputTree;
			return 0;
		}
		if(metadata == this) {
			printWarn << "cannot merge metadata with itself" << endl;
			delete outputTree;
			return 0;
		}
		if(first) {
			first = false;
			if((productionKinematicsParticleNames().empty()) and
			   (decayKinematicsParticleNames().empty()))
			{
				setProductionKinematicsParticleNames(metadata->productionKinematicsParticleNames());
				setDecayKinematicsParticleNames(metadata->decayKinematicsParticleNames());
			}
			if(binningMap().empty()) {
				setBinningMap(metadata->binningMap());
			}
			if(additionalSavedVariableLables().empty()) {
				setAdditionalSavedVariableLables(metadata->additionalSavedVariableLables());
				additionalSavedVariables.resize(additionalSavedVariableLables().size(), 0.);
				for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
					if(inputTree->SetBranchAddress(additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i]) < 0)
					{
						printWarn << "could not set address for branch '" << additionalSavedVariableLables()[i] << "'." << endl;
						delete outputTree;
						return 0;
					}
				}
			}
		}
		appendToUserString(metadata->userString());
		if(productionKinematicsParticleNames() != metadata->productionKinematicsParticleNames()) {
			printWarn << "particle names of production kinematics differ." << endl;
			delete outputTree;
			return 0;
		}
		if(decayKinematicsParticleNames() != metadata->decayKinematicsParticleNames()) {
			printWarn << "particle names of decay kinematics differ." << endl;
			delete outputTree;
			return 0;
		}
		if(binningMap() != metadata->binningMap()) {
			printWarn << "binning maps differ." << endl;
			delete outputTree;
			return 0;
		}
		if(inputTree->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
			delete outputTree;
			return 0;
		}
		if(inputTree->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
			delete outputTree;
			return 0;
		}
		for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
			if(inputTree->SetBranchAddress(additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i]) < 0)
			{
				printWarn << "could not set address for branch '" << additionalSavedVariableLables()[i] << "'." << endl;
				delete outputTree;
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
			outputTree->Fill();
		}

	}
	setContentHash(hashor.hash());
	return outputTree;
}
