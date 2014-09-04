
#ifndef EVENTMETADATA_H
#define EVENTMETADATA_H

#include <map>

#include <TObject.h>

class TTree;


namespace rpwa {

	class eventMetadata : public TObject {
		friend class eventFileWriter;
		friend class eventStorage;

	  private:
		typedef std::pair<double, double> rangePairType;
		typedef std::map<std::string, rangePairType> binningMapType;

	  public:

		eventMetadata();
		~eventMetadata();

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const binningMapType& binningMap() const { return _binningMap; }
		const std::vector<std::string>& productionKinematicsParticleNames() const { return _productionKinematicsParticleNames; }
		const std::vector<std::string>& decayKinematicsParticleNames() const { return _decayKinematicsParticleNames; }
		const std::vector<std::string>& additionalSavedVariableLables() const { return _additionalSavedVariableLabels; }

		std::string hash(TTree* eventTree, const bool& printProgress = false) const;

		std::ostream& print(std::ostream& out) const;

		Long64_t Merge(TCollection* list, Option_t* option = "");     // throws an exception
		TTree* merge(const std::vector<std::pair<const rpwa::eventMetadata*, TTree*> >& inputData,
		             const int& splitlevel = 99,
		             const int& buffsize = 256000);                     // actually works

		static const std::string objectNameInFile;
		static const std::string eventTreeName;
		static const std::string productionKinematicsMomentaBranchName;
		static const std::string decayKinematicsMomentaBranchName;

	  private:

		void setUserString(const std::string& userString) { _userString = userString; }
		void appendToUserString(const std::string& userString,
		                        const std::string& delimiter = ", ");

		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setProductionKinematicsParticleNames(const std::vector<std::string>& productionKinematicsParticleNames);
		void setDecayKinematicsParticleNames(const std::vector<std::string>& decayKinematicsParticleNames);
		void setAdditionalSavedVariableLables(std::vector<std::string> labels) { _additionalSavedVariableLabels = labels; }

		void setBinningVariableLabels(const std::vector<std::string>& labels);
		void setBinningVariableRange(const std::string& label, const rangePairType& range);
		void setBinningMap(const binningMapType& binningMap);

		std::string _userString;
		std::string _contentHash;

		std::vector<std::string> _productionKinematicsParticleNames;
		std::vector<std::string> _decayKinematicsParticleNames;

		binningMapType _binningMap;

		std::vector<std::string> _additionalSavedVariableLabels;

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
