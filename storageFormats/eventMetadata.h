
#ifndef EVENTMETADATA_H
#define EVENTMETADATA_H

#include <map>

#include <TObject.h>

class TTree;
class TFile;


namespace rpwa {

	class eventMetadata : public TObject {
		friend class eventFileWriter;

	  private:
		typedef std::pair<double, double> rangePairType;
		typedef std::map<std::string, rangePairType> binningMapType;

	  public:

		enum eventsTypeEnum {
			OTHER,
			REAL,
			GENERATED,
			ACCEPTED
		};

		~eventMetadata();

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const eventsTypeEnum& eventsType() const { return _eventsType; }
		const binningMapType& binningMap() const { return _binningMap; }
		const std::vector<std::string>& productionKinematicsParticleNames() const { return _productionKinematicsParticleNames; }
		const std::vector<std::string>& decayKinematicsParticleNames() const { return _decayKinematicsParticleNames; }
		const std::vector<std::string>& additionalSavedVariableLables() const { return _additionalSavedVariableLabels; }

		std::string recalculateHash(const bool& printProgress = false) const;

		std::ostream& print(std::ostream& out) const;

		Long64_t Merge(TCollection* list, Option_t* option = "");     // throws an exception
		static eventMetadata* merge(const std::vector<const rpwa::eventMetadata*>& inputData,
		                            const int& splitlevel = 99,
		                            const int& buffsize = 256000);                     // actually works

		TTree* eventTree() const { return _eventTree; } // changing this tree is not allowed (it should be const, but then you can't read it...)

		static const eventMetadata* readEventFile(TFile* inputFile, const bool& quiet = false);

		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) { return ((const eventMetadata*)this)->Write(name, option, bufsize); }
		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

		static const std::string objectNameInFile;
		static const std::string eventTreeName;
		static const std::string productionKinematicsMomentaBranchName;
		static const std::string decayKinematicsMomentaBranchName;

#if defined(__CINT__) || defined(__CLING__) || defined(G__DICTIONARY)
	// root needs a public default constructor
	  public:
#else
	  private:
#endif

		eventMetadata();

	  private:

		void setUserString(const std::string& userString) { _userString = userString; }
		void appendToUserString(const std::string& userString,
		                        const std::string& delimiter = ", ");

		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setEventsType(const eventsTypeEnum& eventsType) { _eventsType = eventsType; }
		void setProductionKinematicsParticleNames(const std::vector<std::string>& productionKinematicsParticleNames);
		void setDecayKinematicsParticleNames(const std::vector<std::string>& decayKinematicsParticleNames);
		void setAdditionalSavedVariableLables(std::vector<std::string> labels) { _additionalSavedVariableLabels = labels; }

		void setBinningVariableLabels(const std::vector<std::string>& labels);
		void setBinningVariableRange(const std::string& label, const rangePairType& range);
		void setBinningMap(const binningMapType& binningMap);

		static std::string getStringForEventsType(const eventsTypeEnum& type);

		std::string _userString;
		std::string _contentHash;
		eventsTypeEnum _eventsType;

		std::vector<std::string> _productionKinematicsParticleNames;
		std::vector<std::string> _decayKinematicsParticleNames;

		binningMapType _binningMap;

		std::vector<std::string> _additionalSavedVariableLabels;

		mutable TTree* _eventTree; //!

		ClassDef(eventMetadata, 1);

	}; // class eventMetadata


	inline
	std::ostream&
	operator <<(std::ostream&          out,
	            const eventMetadata&    metadata)
	{
		return metadata.print(out);
	}

} // namespace rpwa

#endif
