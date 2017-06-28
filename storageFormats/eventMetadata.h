
#ifndef EVENTMETADATA_H
#define EVENTMETADATA_H

#include <TObject.h>

#include "multibinTypes.h"

class TFile;
class TTree;


namespace rpwa {

	class eventMetadata : public TObject {
		friend class eventFileWriter;

	  public:

		enum eventsTypeEnum {
			OTHER,
			REAL,
			GENERATED,
			ACCEPTED
		};

		~eventMetadata();

		bool operator==(const eventMetadata& rhs) const;
		bool operator!=(const eventMetadata& rhs) const { return not (*this == rhs); }

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const eventsTypeEnum& eventsType() const { return _eventsType; }
		const rpwa::multibinBoundariesType& multibinBoundaries() const { return _multibinBoundaries; }
		const std::vector<std::string>& productionKinematicsParticleNames() const { return _productionKinematicsParticleNames; }
		const std::vector<std::string>& decayKinematicsParticleNames() const { return _decayKinematicsParticleNames; }
		const std::vector<std::string>& additionalSavedVariableLables() const { return _additionalSavedVariableLabels; }

		std::string recalculateHash(const bool& printProgress = false) const;

		std::ostream& print(std::ostream& out) const;

		Long64_t Merge(TCollection* list, Option_t* option = "");     // throws an exception
		static eventMetadata* merge(const std::vector<const rpwa::eventMetadata*>& inputData,
		                            const bool mergeDiffMeta = false,
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
		void setBinningVariableRange(const std::string& label, const rpwa::boundaryType& range);
		void setMultibinBoundaries(const rpwa::multibinBoundariesType& multibinBoundaries);

		static std::string getStringForEventsType(const eventsTypeEnum& type);

		std::string _userString;
		std::string _contentHash;
		eventsTypeEnum _eventsType;

		std::vector<std::string> _productionKinematicsParticleNames;
		std::vector<std::string> _decayKinematicsParticleNames;

		rpwa::multibinBoundariesType _multibinBoundaries;

		std::vector<std::string> _additionalSavedVariableLabels;

		mutable TTree* _eventTree; //!

		ClassDef(eventMetadata, 2);

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
