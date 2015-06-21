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
#include "waveDescription_py.h"
#include "amplitudeTreeLeaf_py.h"
#include "ampIntegralMatrix_py.h"
#include "phaseSpaceIntegral_py.h"
#include "nBodyPhaseSpaceGen_py.h"
#include "randomNumberGenerator_py.h"
#include "generatorManager_py.h"
#include "generator_py.h"
#include "generatorParameters_py.h"
#include "generatorPickerFunctions_py.h"
#include "beamAndVertexGenerator_py.h"
#include "complexMatrix_py.h"
#include "fitResult_py.h"
#include "physUtils_py.h"
#include "reportingUtilsEnvironment_py.h"
#include "eventFileWriter_py.h"
#include "eventMetadata_py.h"
#include "amplitudeFileWriter_py.h"
#include "amplitudeMetadata_py.h"
#include "calcAmplitude_py.h"
#include "pwaFit_py.h"
#include "pwaLikelihood_py.h"
#ifdef USE_NLOPT
#include "pwaNloptFit_py.h"
#endif

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
	rpwa::py::exportWaveDescription();
	rpwa::py::exportAmplitudeTreeLeaf();
	rpwa::py::exportAmpIntegralMatrix();
	rpwa::py::exportPhaseSpaceIntegral();
	rpwa::py::exportNBodyPhaseSpaceGen();
	rpwa::py::exportRandomNumberGenerator();
	rpwa::py::exportGeneratorManager();
	rpwa::py::exportGenerator();
	rpwa::py::exportGeneratorParameters();
	rpwa::py::exportGeneratorPickerFunctions();
	rpwa::py::exportBeamAndVertexGenerator();
	rpwa::py::exportComplexMatrix();
	rpwa::py::exportFitResult();
	rpwa::py::exportPhysUtils();
	rpwa::py::exportReportingUtilsEnvironment();
	rpwa::py::exportEventFileWriter();
	rpwa::py::exportEventMetadata();
	rpwa::py::exportAmplitudeFileWriter();
	rpwa::py::exportAmplitudeMetadata();
	rpwa::py::exportCalcAmplitude();
	rpwa::py::exportPwaLikelihood();
	rpwa::py::exportPwaFit();
#ifdef USE_NLOPT
	rpwa::py::exportPwaNloptFit();
#endif

}
