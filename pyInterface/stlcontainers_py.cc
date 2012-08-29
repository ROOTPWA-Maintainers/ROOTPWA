#include<stlcontainers_py.h>

#include<particleProperties.h>

namespace bp = boost::python;

void rpwa::py::exportStdPairs() {

	typedef std::pair< const std::string, rpwa::particleProperties > stdpair_int_particleProperties;
	bp::class_< stdpair_int_particleProperties >("__stdpair_int_particleProperties")
		.add_property("first", &stdpair_int_particleProperties::first)
		.def_readwrite("second", &stdpair_int_particleProperties::second);

};


std::set<std::string> rpwa::py::converBPObjectToStrSet(boost::python::object list) {

			bp::list pyList = bp::extract<bp::list>(list);
			std::set<std::string> set;
			for(unsigned int i = 0; i < bp::len(pyList); ++i) {
				std::string entry = bp::extract<std::string>(pyList[i]);
				set.insert(entry);
			}
			return set;

};
