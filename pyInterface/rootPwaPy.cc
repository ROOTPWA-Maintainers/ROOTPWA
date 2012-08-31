#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "stlContainers_py.h"
#include "rootConverters_py.h"
#include "particleProperties_py.h"
#include "particleDataTable_py.h"
#include "particle_py.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(libRootPwaPy){

	rpwa::py::exportStdPairs();
	rpwa::py::exportRootConverters();
	rpwa::py::exportParticlePropertiesVector();
	rpwa::py::exportParticleProperties();
	rpwa::py::exportParticleDataTable();
	rpwa::py::exportParticle();

}


