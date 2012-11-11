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

		void exportStlContainers();

		std::set     <std::string> convertBPObjectToStrSet     (boost::python::object list);
		std::multiset<std::string> convertBPObjectToStrMultiSet(boost::python::object list);

	}
}

#endif
