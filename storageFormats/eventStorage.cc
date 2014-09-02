
#include "eventStorage.h"

#include <boost/progress.hpp>

#include <TClonesArray.h>

#include "reportingUtils.hpp"
#include "hashCalculator.h"

using namespace std;
using namespace rpwa;


const std::string rpwa::eventStorage::objectNameInFile = "eventStorageROOTPWA";
const std::string rpwa::eventStorage::eventTreeName = "rootPwaEvtTree";
const std::string rpwa::eventStorage::productionKinematicsMomentaBranchName = "prodKinMomenta";
const std::string rpwa::eventStorage::decayKinematicsMomentaBranchName = "decayKinMomenta";


rpwa::eventStorage::eventStorage()
	: storage()
{

}


rpwa::eventStorage::~eventStorage()
{

}


bool rpwa::eventStorage::setBranchAddressProductionKinematicsMomenta(TClonesArray*& productionKinematicsMomenta,
                                                                     const int&     splitlevel,
                                                                     const int&     buffsize)
{
	if(data()->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
		printWarn << "could not set address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
		return false;
	}
	return true;
}



bool rpwa::eventStorage::setBranchAddressDecayKinematicsMomenta(TClonesArray*& decayKinematicsMomenta,
                                                                const int&     splitlevel,
                                                                const int&     buffsize)
{
	if(data()->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
		printWarn << "could not set address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
		return false;
	}
	return true;
}


bool rpwa::eventStorage::setBranchAddressAdditionalVariables(std::vector<double>& values,
                                                             const int&           splitlevel,
                                                             const int&           buffsize)
{
	if(values.size() != _metadata.additionalSavedVariableLables().size()) {
		printWarn << "size mismatch between argument vector (" << values.size()
		          << ") and stored vector (" << _metadata.additionalSavedVariableLables().size()
		          << ")." << endl;
		return false;
	}
	for(unsigned int i = 0; i < values.size(); ++i) {
		if(data()->SetBranchAddress(_metadata.additionalSavedVariableLables()[i].c_str(), &(values[i])) < 0) {
			printWarn << "could not set address for branch '" << _metadata.additionalSavedVariableLables()[i] << "'." << endl;
			return false;
		}
	}
	return true;
}


string rpwa::eventStorage::hash(const bool& printProgress)
{
	return startHash(printProgress).hash();
}


hashCalculator rpwa::eventStorage::startHash(const bool& printProgress)
{
	TClonesArray* productionKinematicsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	hashCalculator hashor;
	if(not setBranchAddressProductionKinematicsMomenta(productionKinematicsMomenta))
	{
		return hashor;
	}
	if(not setBranchAddressDecayKinematicsMomenta(decayKinematicsMomenta)) {
		return hashor;
	}
	vector<double> additionalVariables(_metadata.additionalSavedVariableLables().size(), 0.);
	for(unsigned int i = 0; i < additionalVariables.size(); ++i) {
		data()->SetBranchAddress(_metadata.additionalSavedVariableLables()[i].c_str(), &additionalVariables[i]);
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(data()->GetEntries(), cout, "") : 0;
	for(long eventNumber = 0; eventNumber < data()->GetEntries(); ++eventNumber) {
		data()->GetEntry(eventNumber);
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
	return hashor;
}


Long64_t rpwa::eventStorage::Merge(TCollection* list, Option_t* option)
{
	TIter next(list);
	eventStorage* storage = 0;
	bool first = true;
	hashCalculator hashor;
	if(data()->GetEntries() > 0) {
		hashor = startHash();
	}
	TClonesArray* productionKinematicsMomenta = 0;
	TClonesArray* decayKinematicsMomenta = 0;
	if(data()->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
		printWarn << "could not set branch address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
		return -1;
	}
	if(data()->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
			return -1;
	}
	vector<double> additionalSavedVariables;
	while((storage = dynamic_cast<eventStorage*>(next()))) {
		if(not storage) {
			printWarn << "input object of the wrong type when merging." << endl;
			return -1;
		}
		if(storage == this) {
			continue;
		}
		if(first) {
			first = false;
			if((_metadata.productionKinematicsParticleNames().empty()) and
			   (_metadata.decayKinematicsParticleNames().empty()))
			{
				_metadata.setProductionKinematicsParticleNames(storage->metadata().productionKinematicsParticleNames());
				_metadata.setDecayKinematicsParticleNames(storage->metadata().decayKinematicsParticleNames());
			}
			if(_metadata.binningMap().empty()) {
				_metadata.setBinningMap(storage->metadata().binningMap());
			}
			if(_metadata.additionalSavedVariableLables().empty()) {
				_metadata.setAdditionalSavedVariableLables(storage->metadata().additionalSavedVariableLables());
				additionalSavedVariables.resize(_metadata.additionalSavedVariableLables().size(), 0.);
				for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
					if(data()->SetBranchAddress(_metadata.additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i]) < 0)
					{
						printWarn << "could not set address for branch '" << _metadata.additionalSavedVariableLables()[i] << "'." << endl;
						return -1;
					}
				}
			}
		}
		_metadata.appendToUserString(storage->metadata().userString());
		if(_metadata.productionKinematicsParticleNames() != storage->metadata().productionKinematicsParticleNames()) {
			printWarn << "particle names of production kinematics differ." << endl;
			return -1;
		}
		if(_metadata.decayKinematicsParticleNames() != storage->metadata().decayKinematicsParticleNames()) {
			printWarn << "particle names of decay kinematics differ." << endl;
			return -1;
		}
		if(_metadata.binningMap() != storage->metadata().binningMap()) {
			printWarn << "binning maps differ." << endl;
			return -1;
		}
		if(storage->data()->SetBranchAddress(productionKinematicsMomentaBranchName.c_str(), &productionKinematicsMomenta) < 0) {
			printWarn << "could not set branch address for branch '" << productionKinematicsMomentaBranchName << "'." << endl;
			return -1;
		}
		if(storage->data()->SetBranchAddress(decayKinematicsMomentaBranchName.c_str(), &decayKinematicsMomenta) < 0) {
				printWarn << "could not set branch address for branch '" << decayKinematicsMomentaBranchName << "'." << endl;
				return -1;
		}
		for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
			if(storage->data()->SetBranchAddress(_metadata.additionalSavedVariableLables()[i].c_str(), &additionalSavedVariables[i]) < 0)
			{
				printWarn << "could not set address for branch '" << _metadata.additionalSavedVariableLables()[i] << "'." << endl;
				return -1;
			}
		}
		for(long eventNumber = 0; eventNumber < storage->data()->GetEntries(); ++eventNumber) {
			storage->data()->GetEntry(eventNumber);
			for(int i = 0; i < productionKinematicsMomenta->GetEntries(); ++i) {
				hashor.Update(*((TVector3*)(*productionKinematicsMomenta)[i]));
			}
			for(int i = 0; i < decayKinematicsMomenta->GetEntries(); ++i) {
				hashor.Update(*((TVector3*)(*decayKinematicsMomenta)[i]));
			}
			for(unsigned int i = 0; i < additionalSavedVariables.size(); ++i) {
				hashor.Update(additionalSavedVariables[i]);
			}
			data()->Fill();
		}

	}
	_metadata.setContentHash(hashor.hash());
	return data()->GetEntries();
}
