
#ifndef EVENTMETADATA_H
#define EVENTMETADATA_H

#include <map>

#include <TObject.h>


namespace rpwa {

	class eventMetadata : public TObject {
		friend class eventFileWriter;

	  private:
		typedef std::pair<double, double> rangePairType;
		typedef std::map<std::string, rangePairType> binningMapType;

	  public:

		eventMetadata();
		~eventMetadata();

		const std::string& userString() const { return _userString; }
		const std::string& contentHash() const { return _contentHash; }
		const binningMapType& getBinningMap() const { return _binningMap; }
		const std::vector<std::string>& initalStateParticleNames() const { return _initialStateParticleNames; }
		const std::vector<std::string>& finalStateParticleNames() const { return _finalStateParticleNames; }
		const std::vector<std::string>& additionalSavedVariableLables() const { return _additionalSavedVariableLabels; }

		std::ostream& print(std::ostream& out) const;


	  private:

		void setUserString(const std::string& userString) { _userString = userString; }
		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setInitialStateParticleNames(const std::vector<std::string>& initialStateParticleNames);
		void setFinalStateParticleNames(const std::vector<std::string>& finalStateParticleNames);
		void setAdditionalSavedVariableLables(std::vector<std::string> labels) { _additionalSavedVariableLabels = labels; }

		void setBinningVariableLabels(const std::vector<std::string>& labels);
		void setBinningVariableRange(const std::string& label, const rangePairType& range);
		void setBinningMap(const binningMapType& binningMap);

		std::string _userString;
		std::string _contentHash;

		std::vector<std::string> _initialStateParticleNames;
		std::vector<std::string> _finalStateParticleNames;

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
