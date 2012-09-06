#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "stlContainers_py.h"
#include "rootConverters_py.h"
#include "particleProperties_py.h"
#include "particleDataTable_py.h"
#include "particle_py.h"
#include "interactionVertex_py.h"
#include "fsVertex_py.h"
#include "massDependence_py.h"
#include "isobarDecayVertex_py.h"
#include "productionVertex_py.h"
#include "diffractiveDissVertex_py.h"
#include "decayTopology_py.h"
#include "isobarDecayTopology_py.h"
#include "isobarAmplitude_py.h"
#include "isobarCanonicalAmplitude_py.h"
#include "isobarHelicityAmplitude_py.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(libRootPwaPy){

	rpwa::py::exportStlContainers();
	rpwa::py::exportParticleProperties();
	rpwa::py::exportParticleDataTable();
	rpwa::py::exportParticle();
	rpwa::py::exportInteractionVertex();
	rpwa::py::exportFsVertex();
	rpwa::py::exportMassDependence();
	rpwa::py::exportIsobarDecayVertex();
	rpwa::py::exportProductionVertex();
	rpwa::py::exportDiffractiveDissVertex();
	rpwa::py::exportDecayTopology();
	rpwa::py::exportIsobarDecayTopology();
	rpwa::py::exportIsobarAmplitude();
	rpwa::py::exportIsobarCanonicalAmplitude();
	rpwa::py::exportIsobarHelicityAmplitude();

}

