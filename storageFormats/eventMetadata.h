
#ifndef EVENTMETADATA_H
#define EVENTMETADATA_H

#include <TObject.h>

#include "multibinTypes.h"

class TFile;
class TTree;


namespace rpwa {

	class additionalTreeVariables;
	class hashCalculator;

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

		/***
		 * Compare the content of both metadata. Auxiliary information (auxString, auxValues) are not compared.
		 */
		bool operator==(const eventMetadata& rhs) const;
		bool operator!=(const eventMetadata& rhs) const { return not (*this == rhs); }

		const std::string& auxString() const { return _auxString; }
		const std::string& contentHash() const { return _contentHash; }
		const eventsTypeEnum& eventsType() const { return _eventsType; }
		const rpwa::multibinBoundariesType& multibinBoundaries() const { return _multibinBoundaries; }
		const std::vector<std::string>& productionKinematicsParticleNames() const { return _productionKinematicsParticleNames; }
		const std::vector<std::string>& decayKinematicsParticleNames() const { return _decayKinematicsParticleNames; }
		const std::vector<std::string>& additionalTreeVariableNames() const { return _additionalTreeVariableNames; }
		const std::map<std::string, double>& auxValues() const { return _auxValues; }
		double auxValue(const std::string& name) const { return _auxValues.at(name); }
		void setAuxValue(const std::string& name, const double value) { _auxValues[name] = value; }
		bool hasAuxValue(const std::string& name) const { return _auxValues.find(name) != _auxValues.end(); }

		/***
		 * Update the given hashor with the default information used to build the hash of the meta data
		 * \return true if update was successful
		 */
		bool updateHashor(hashCalculator& hashor, const bool& printProgress = false) const;
		/***
		 * Recalculate the hash for this object. Auxiliary information (auxString, auxValues) are not considered.
		 */
		std::string recalculateHash(const bool& printProgress = false) const;

		std::ostream& print(std::ostream& out) const;

		Long64_t Merge(TCollection* list, Option_t* option = "");     // throws an exception
		/***
		 * \par mergeAuxValues Also merge the auxiliary values if they are the same
		 * @return a new eventMetadata, which is the merge of all inputData.
		 *         By default, auxiliary information (auxString and auxValues) are not merged.
		 */
		static eventMetadata* merge(const std::vector<const rpwa::eventMetadata*>& inputData,
		                            const bool mergeBinBoundaries = false,
		                            const bool mergeAuxString = false,
		                            const bool mergeAuxValues = false,
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

		void setAuxString(const std::string& auxString) { _auxString = auxString; }
		void appendToAuxString(const std::string& auxString,
		                       const std::string& delimiter = ", ");

		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setEventsType(const eventsTypeEnum& eventsType) { _eventsType = eventsType; }
		void setProductionKinematicsParticleNames(const std::vector<std::string>& productionKinematicsParticleNames) { _productionKinematicsParticleNames = productionKinematicsParticleNames; }
		void setDecayKinematicsParticleNames(const std::vector<std::string>& decayKinematicsParticleNames) { _decayKinematicsParticleNames = decayKinematicsParticleNames; }
		void setAdditionalTreeVariableNames(const std::vector<std::string>& labels) { _additionalTreeVariableNames = labels; }

		void setBinningVariableLabels(const std::vector<std::string>& labels);
		void setBinningVariableRange(const std::string& label, const rpwa::boundaryType& range) { _multibinBoundaries[label] = range; }
		void setMultibinBoundaries(const rpwa::multibinBoundariesType& multibinBoundaries) { _multibinBoundaries = multibinBoundaries; }

		static std::string getStringForEventsType(const eventsTypeEnum& type);

		std::string _auxString; // the content of this variable is by default not included in the '==' comparison, hash calculation, or merging
		std::string _contentHash;
		eventsTypeEnum _eventsType;

		std::vector<std::string> _productionKinematicsParticleNames;
		std::vector<std::string> _decayKinematicsParticleNames;

		rpwa::multibinBoundariesType _multibinBoundaries;

		std::vector<std::string> _additionalTreeVariableNames;

		std::map<std::string, double> _auxValues; // the content of this variable is by default not included in the '==' comparison, hash calculation, or merging

		mutable TTree* _eventTree; //!

		ClassDef(eventMetadata, 4);

	}; // class eventMetadata


	inline
	std::ostream&
	operator <<(std::ostream&          out,
	            const eventMetadata&    metadata)
	{
		return metadata.print(out);
	}


	/**
	 * \brief handler for additional kinematic variables of event data
	 *
	 * The class represents a set of kinematic variables used for on-the-fly binning,
	 * which are stored as additional tree variables together with the event data.
	 * It provides an interface to connect kinematic variables to tree branches
	 * and to check whether an event is within a kinematic bin.
	 */
	class additionalTreeVariables {

	  public:

		additionalTreeVariables() {}
		additionalTreeVariables(const additionalTreeVariables&) = delete;
		additionalTreeVariables& operator= (const additionalTreeVariables&) = delete;

		double operator[](const std::string& label) const { return _additionalTreeVariables.at(label); }

		/**
		 * reads list of tree variables from provided metaData and connects them to tree branches
		 * \param metaData eventMetadata object, whose additional tree variables branches will be set to this object
		 * \return true if the setting of the branch addresses was successful
		 */
		bool setBranchAddresses(const eventMetadata& metaData);

		/**
		 * checks whether the current additional variables are within the boundaries
		 * x_i in [lower_i, upper_i) for all variables i in the given boundaries
		 */
		bool inBoundaries(const multibinBoundariesType& boundaries) const;

	  private:

		std::map<std::string, double> _additionalTreeVariables;

	};

} // namespace rpwa

#endif
