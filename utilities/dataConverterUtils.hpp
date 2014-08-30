
#ifndef DATACONVERTERUTILS_HPP
#define DATACONVERTERUTILS_HPP

#include <map>

#include <TClonesArray.h>
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <TVector3.h>

#include "reportingUtils.hpp"


namespace rpwa {

	class dataMetadata : public TObject {
		friend class rootpwaDataFileWriter;


	  public:

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const std::map<std::string, std::pair<double, double> >& getBinningMap() const { return _binningMap; }
		const std::vector<std::string>& initalStateParticleNames() const { return _initialStateParticleNames; }
		const std::vector<std::string>& finalStateParticleNames() const { return _finalStateParticleNames; }

		std::ostream& print(std::ostream& out) const
		{
			out << "dataMetadata: " << std::endl
			    << "    userString ...................... '" << _userString << "'"         << std::endl
			    << "    contentHash ..................... '" << _contentHash << "'"        << std::endl
			    << "    initial state particle names: ... "  << _initialStateParticleNames << std::endl
		        << "    final state particle names: ..... "  << _finalStateParticleNames   << std::endl
			    << "    binning map: "                       << std::endl;
			    for(std::map<std::string, std::pair<double, double> >::const_iterator it = _binningMap.begin(); it != _binningMap.end(); ++it) {
					out << "        variable '" << it->first << "' range " << it->second << std::endl;
			    }
			return out;
		}


	  private:

		void setUserString(const std::string& userString) { _userString = userString; }
		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setInitialStateParticleNames(std::vector<std::string> initialStateParticleNames) { _initialStateParticleNames = initialStateParticleNames; }
		void setFinalStateParticleNames(std::vector<std::string> finalStateParticleNames) { _finalStateParticleNames = finalStateParticleNames; }

		void setBinningVariableLabels(std::vector<std::string> labels)
		{
			for(unsigned int i = 0; i < labels.size(); ++i) {
				_binningMap[labels[i]] = std::pair<double, double>(0., 0.);
			}
		}

		void setBinningVariableRange(const std::string& label, const std::pair<double, double>& range)
		{
			_binningMap[label] = range;
		}

		void setBinningMap(const std::map<std::string, std::pair<double, double> > binningMap)
		{
			_binningMap = binningMap;
		}

		std::string _userString;
		std::string _contentHash;

		std::vector<std::string> _initialStateParticleNames;
		std::vector<std::string> _finalStateParticleNames;

		std::map<std::string, std::pair<double, double> > _binningMap;

	}; // class dataMetadata


	inline
	std::ostream&
	operator <<(std::ostream&          out,
	            const dataMetadata&    metadata)
	{
		return metadata.print(out);
	}


	class rootpwaDataFileWriter {

	  public:

		rootpwaDataFileWriter()
			: _initialized(false),
			  _outfile(0),
			  _eventTree(0),
			  _initialStateMomenta(0),
			  _finalStateMomenta(0),
			  _additionalVariablesToSave(),
			  _metadata(),
			  _nmbInitialStateParticles(0),
			  _nmbFinalStateParticles(0) { }

		~rootpwaDataFileWriter()
		{
			reset();
		}

		bool initialize(TFile&                                                   outputFile,                // output file to write the data to (user keeps ownership!)
		                const std::string&                                       userString,                // some arbitrary string to identify this data file
		                const std::vector<std::string>&                          initialStateParticleNames, // particle names of initial state particles (has to be the same order as the particles appear in the data!)
		                const std::vector<std::string>&                          finalStateParticleNames,   // particle names of final state particles (has to be the same order as the particles appear in the data!)
		                const std::map<std::string, std::pair<double, double> >& binningMap,                // binning variable map with content "label" -> (lowerBound, upperBound) describing which bin these data belong to
		                const std::vector<std::string>&                          additionalVariableLabels,  // Labels for any additional information which is stored (as double) and can later be used for binning
		                const std::string&              eventTreeName = "rootPwaEvtTree", // name for the event tree in the output file
		                const std::string&              initialStateMomentaBranchName = "prodKinMomenta",  // branch name where the initial state particles are stored
		                const std::string&              finalStateMomentaBranchName   = "decayKinMomenta", // branch name where the final state particles are stored
		                const int&                      splitlevel = 99,
		                const int&                      buffsize = 256000)
		{
			if(_initialized) {
				printWarn << "trying to initialize when already initialized" << std::endl;
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
			_additionalVariablesToSave = std::vector<double>(additionalVariableLabels.size(), 0.);
			for(unsigned int i = 0; i < additionalVariableLabels.size(); ++i) {
				std::stringstream strStr;
				strStr << additionalVariableLabels[i] << "/D";
				_eventTree->Branch(additionalVariableLabels[i].c_str(), &_additionalVariablesToSave[i], strStr.str().c_str());
			}

			_initialized = true;
			return _initialized;
		}

		void addEvent(const std::vector<TVector3>& initialStateMomenta,
		              const std::vector<TVector3>& finalStateMomenta,
		              const std::vector<double>& additionalVariablesToSave)
		{
			// TODO: calculate hash
			if(initialStateMomenta.size() != _nmbInitialStateParticles) {
				printErr << "received unexpected number of initial state particles (got "
				         << initialStateMomenta.size() << ", expected "
				         << _nmbInitialStateParticles << "). Aborting..." << std::endl;
				throw;
			}
			if(finalStateMomenta.size() != _nmbFinalStateParticles) {
				printErr << "received unexpected number of final state particles (got "
				         << finalStateMomenta.size() << ", expected "
				         << _nmbFinalStateParticles << "). Aborting..." << std::endl;
				throw;
			}
			if(additionalVariablesToSave.size() != _additionalVariablesToSave.size()) {
				printErr << "recieved unexpected number of additional variables (got "
				         << additionalVariablesToSave.size() << ", expected "
				         << _additionalVariablesToSave.size() << "). Aborting..." << std::endl;
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

		bool finalize() {
			if(not _initialized) {
				printWarn << "trying to finalize when not initialized" << std::endl;
				return false;
			}
			_outfile->cd();
			_metadata.Write();
			_outfile->Write();
			_outfile->Close();
			reset();
			return true;
		}

		void reset() {
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

		const bool& initialized() { return _initialized; }

	  private:

		bool _initialized;
		TFile* _outfile;
		TTree* _eventTree;
		TClonesArray* _initialStateMomenta;
		TClonesArray* _finalStateMomenta;
		std::vector<double> _additionalVariablesToSave;
		dataMetadata _metadata;
		unsigned int _nmbInitialStateParticles;
		unsigned int _nmbFinalStateParticles;

	}; // rootpwaDataFileWriter

} // namespace rpwa

#endif
