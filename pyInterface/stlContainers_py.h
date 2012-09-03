#ifndef STLCONTAINERS_PY_H
#define STLCONTAINERS_PY_H

#include "boost/python.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include<map>
#include<set>
#include<utility>
#include<vector>

namespace rpwa {
	namespace py {

		void exportStdPairs();

		void exportParticlePropertiesVector();
		void exportParticlePtrVector();

		std::set<std::string> converBPObjectToStrSet(boost::python::object list);

	}
}

#endif
