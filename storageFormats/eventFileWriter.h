
#ifndef DATACONVERTERUTILS_HPP
#define DATACONVERTERUTILS_HPP

#include <map>

#include "eventMetadata.h"
#include "hashCalculator.h"

class TFile;
class TVector3;
class TTree;

namespace rpwa {

	class eventFileWriter {

	  public:

		eventFileWriter();
		~eventFileWriter();

		bool initialize(TFile&                                                   outputFile,                        // output file to write the data to (user keeps ownership!)
		                const std::string&                                       userString,                        // some arbitrary string to identify this data file
		                const eventMetadata::eventsTypeEnum&                     eventsType,                        // type of events
		                const std::vector<std::string>&                          productionKinematicsParticleNames, // particle names of initial state particles (has to be the same order as the particles appear in the data!)
		                const std::vector<std::string>&                          decayKinematicsParticleNames,      // particle names of final state particles (has to be the same order as the particles appear in the data!)
		                const std::map<std::string, std::pair<double, double> >& binningMap,                        // binning variable map with content "label" -> (lowerBound, upperBound) describing which bin these data belong to
		                const std::vector<std::string>&                          additionalVariableLabels,          // Labels for any additional information which is stored (as double) and can later be used for binning
		                const int&                                               splitlevel = 99,
		                const int&                                               buffsize = 256000);

		void addEvent(const std::vector<TVector3>& productionKinematicsMomenta,
		              const std::vector<TVector3>& decayKinematicsMomenta,
		              const std::vector<double>&   additionalVariablesToSave = std::vector<double>());

		void addEvent(const TClonesArray&        productionKinematicsMomenta,
		              const TClonesArray&        decayKinematicsMomenta,
		              const std::vector<double>& additionalVariablesToSave = std::vector<double>());


		bool finalize();

		void reset();

		const bool& initialized() { return _initialized; }

	  private:

		bool _initialized;
		TFile* _outputFile;
		eventMetadata _metadata;
		TClonesArray* _productionKinematicsMomenta;
		TClonesArray* _decayKinematicsMomenta;
		std::vector<double> _additionalVariablesToSave;
		unsigned int _nmbProductionKinematicsParticles;
		unsigned int _nmbDecayKinematicsParticles;
		hashCalculator _hashCalculator;

	}; // rootpwaDataFileWriter

} // namespace rpwa

#endif
