
#ifndef DATACONVERTERUTILS_HPP
#define DATACONVERTERUTILS_HPP

#include <map>

#include <TMD5.h>

#include "eventMetadata.h"

class TFile;
class TVector3;
class TTree;

namespace rpwa {

	class md5Wrapper : public TMD5 {

	  public:

		md5Wrapper()
			: TMD5() { }

		void Update(const double& value);
		void Update(const TVector3& vector);

		std::string hash() {
			TMD5::Final();
			return TMD5::AsString();
		}

	};

	class eventFileWriter {

	  public:

		eventFileWriter();
		~eventFileWriter();

		bool initialize(TFile&                                                   outputFile,                // output file to write the data to (user keeps ownership!)
		                const std::string&                                       userString,                // some arbitrary string to identify this data file
		                const std::vector<std::string>&                          initialStateParticleNames, // particle names of initial state particles (has to be the same order as the particles appear in the data!)
		                const std::vector<std::string>&                          finalStateParticleNames,   // particle names of final state particles (has to be the same order as the particles appear in the data!)
		                const std::map<std::string, std::pair<double, double> >& binningMap,                // binning variable map with content "label" -> (lowerBound, upperBound) describing which bin these data belong to
		                const std::vector<std::string>&                          additionalVariableLabels,  // Labels for any additional information which is stored (as double) and can later be used for binning
		                const std::string&              eventTreeName = "rootPwaEvtTree", // name for the event tree in the output file
		                const std::string&              initialStateMomentaBranchName = "prodKinMomenta",  // branch name where the initial state particles are stored
		                const std::string&              finalStateMomentaBranchName   = "decayKinMomenta", // branch name where the final state particles are stored
		                const std::string&              metadataName = "dataMetadata",                     // name under which the metadata object is saved
		                const int&                      splitlevel = 99,
		                const int&                      buffsize = 256000);

		void addEvent(const std::vector<TVector3>& initialStateMomenta,
		              const std::vector<TVector3>& finalStateMomenta,
		              const std::vector<double>& additionalVariablesToSave = std::vector<double>());


		bool finalize();

		void reset();

		const bool& initialized() { return _initialized; }

		static std::string calculateHash(TTree* eventTree,
		                                 const std::vector<std::string> additionalVariableLabels      = std::vector<std::string>(),
		                                 const bool&                    printProgress                 = false,
		                                 const std::string&             initialStateMomentaBranchName = "prodKinMomenta",
		                                 const std::string&             finalStateMomentaBranchName   = "decayKinMomenta");

	  private:

		bool _initialized;
		TFile* _outfile;
		TTree* _eventTree;
		TClonesArray* _initialStateMomenta;
		TClonesArray* _finalStateMomenta;
		std::vector<double> _additionalVariablesToSave;
		eventMetadata _metadata;
		unsigned int _nmbInitialStateParticles;
		unsigned int _nmbFinalStateParticles;
		md5Wrapper _md5Calculator;

	}; // rootpwaDataFileWriter

} // namespace rpwa

#endif
