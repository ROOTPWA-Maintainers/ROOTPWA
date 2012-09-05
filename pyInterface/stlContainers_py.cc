#include "stlContainers_py.h"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include<particle.h>
#include<particleProperties.h>

namespace bp = boost::python;

void rpwa::py::exportStlContainers() {

	// std::pair<std::string, rpwa::particleProperties>
	typedef std::pair< const std::string, rpwa::particleProperties > stdpair_int_particleProperties;
	bp::class_< stdpair_int_particleProperties >("__stdpair_int_particleProperties")
		.add_property("first", &stdpair_int_particleProperties::first)
		.def_readwrite("second", &stdpair_int_particleProperties::second);

	// std::vector<rpwa::particleProperties>
	bp::class_<std::vector<rpwa::particleProperties> >("__vector_particleProperties")
		.def(
			bp::vector_indexing_suite<std::vector<rpwa::particleProperties> >()
		);

	// std::vector<rpwa::particlePtr>
	bp::class_<std::vector<rpwa::particlePtr> >("__vector_particlePtr")
		.def(
			bp::vector_indexing_suite<std::vector<rpwa::particlePtr> >()
		);
};

std::set<std::string> rpwa::py::converBPObjectToStrSet(bp::object list) {

			bp::list pyList = bp::extract<bp::list>(list);
			std::set<std::string> set;
			for(unsigned int i = 0; i < bp::len(pyList); ++i) {
				std::string entry = bp::extract<std::string>(pyList[i]);
				set.insert(entry);
			}
			return set;

};
