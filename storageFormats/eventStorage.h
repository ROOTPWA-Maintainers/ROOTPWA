
#ifndef RPWA_EVENTSTORAGE_H
#define RPWA_EVENTSTORAGE_H

#include "storage.h"
#include "eventMetadata.h"
#include "hashCalculator.h"

namespace rpwa {

	class eventStorage : public storage<eventMetadata>
	{

	  public:

		eventStorage();
		virtual ~eventStorage();

		bool setBranchAddressProductionKinematicsMomenta(TClonesArray*& productionKinematicsMomenta,
		                                                 const int&     splitlevel = 99,
		                                                 const int&     buffsize   = 256000);
		bool setBranchAddressDecayKinematicsMomenta(TClonesArray*& decayKinematicsMomenta,
		                                            const int&     splitlevel = 99,
		                                            const int&     buffsize   = 256000);
		bool setBranchAddressAdditionalVariables(std::vector<double>& values,
		                                         const int&           splitlevel = 99,
		                                         const int&           buffsize   = 256000);

		std::string hash(const bool& printProgress = false);

		Long64_t Merge(TCollection* list, Option_t* option = "");

		static const std::string objectNameInFile;
		static const std::string eventTreeName;
		static const std::string productionKinematicsMomentaBranchName;
		static const std::string decayKinematicsMomentaBranchName;

	  private:

		hashCalculator startHash(const bool& printProgress = false);

		ClassDef(eventStorage, 1);

	};

} // namespace rpwa

#endif
