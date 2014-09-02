
#ifndef RPWA_EVENTSTORAGE_H
#define RPWA_EVENTSTORAGE_H

#include "storage.h"
#include "eventMetadata.h"

namespace rpwa {

	class eventStorage : public storage<eventMetadata>
	{

		Long64_t Merge(TCollection* list, Option_t* option = "");

	};

} // namespace rpwa

#endif
