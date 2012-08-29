#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "particleProperties_py.h"
#include "particleDataTable_py.h"
#include "stlcontainers_py.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(libRootPwaPy){

	rpwa::py::exportStdPairs();
	rpwa::py::exportParticleProperties();
	rpwa::py::exportParticleDataTable();

}


