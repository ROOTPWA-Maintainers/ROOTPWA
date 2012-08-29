#ifndef PARTICLEDATATABLE_PY_H
#define PARTICLEDATATABLE_PY_H

#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "particleDataTable.h"

namespace rpwa {
	namespace py {
		void exportParticleDataTable();
	}
}

#endif
