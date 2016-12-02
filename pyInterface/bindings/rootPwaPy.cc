#include <boost/python.hpp>

#include <RVersion.h>
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 6, 0)
#include <TThread.h>
#else
#include <TROOT.h>
#endif

// set up numpy for usage in this module
#define PY_ARRAY_UNIQUE_SYMBOL RPWA_PyArray_API
#include <numpy/arrayobject.h>

// pyUtils
#include "rootConverters_py.h"
#include "stlContainers_py.h"

// decayAmplitude
#include "ampIntegralMatrix_py.h"
#include "ampIntegralMatrixMetadata_py.h"
#include "decayTopology_py.h"
#include "diffractiveDissVertex_py.h"
#include "fsVertex_py.h"
#include "interactionVertex_py.h"
#include "isobarAmplitude_py.h"
#include "isobarCanonicalAmplitude_py.h"
#include "isobarDecayTopology_py.h"
#include "isobarDecayVertex_py.h"
#include "isobarHelicityAmplitude_py.h"
#include "massDependence_py.h"
#include "phaseSpaceIntegral_py.h"
#include "productionVertex_py.h"
#include "waveDescription_py.h"

// generators
#include "beamAndVertexGenerator_py.h"
#include "generator_py.h"
#include "generatorManager_py.h"
#include "generatorParameters_py.h"
#include "generatorPickerFunctions_py.h"
#ifdef USE_BAT
#include "importanceSampler_py.h"
#endif
#include "modelIntensity_py.h"

// highLevelInterface
#include "calcAmplitude_py.h"
#include "getMassShapes_py.h"
#include "pwaFit_py.h"
#ifdef USE_NLOPT
#include "pwaNloptFit_py.h"
#endif

// nBodyPhaseSpace
#include "nBodyPhaseSpaceGenerator_py.h"
#include "nBodyPhaseSpaceKinematics_py.h"
#include "randomNumberGenerator_py.h"

// partialWaveFit
#include "complexMatrix_py.h"
#include "fitResult_py.h"
#include "partialWaveFitHelper_py.h"
#include "pwaLikelihood_py.h"

// particleData
#include "particle_py.h"
#include "particleDataTable_py.h"
#include "particleProperties_py.h"

// storageFormats
#include "amplitudeFileWriter_py.h"
#include "amplitudeMetadata_py.h"
#include "amplitudeTreeLeaf_py.h"
#include "eventFileWriter_py.h"
#include "eventMetadata_py.h"
#include "hashCalculator_py.h"

// utilities
#include "physUtils_py.h"
#include "reportingUtilsEnvironment_py.h"


BOOST_PYTHON_MODULE(libRootPwaPy){

	// enable multithreading of Python so we can release the GIL whereever
	// appropriate
	PyEval_InitThreads();

	// also enable multithreading of ROOT
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 6, 0)
	TThread::Initialize();
#else
	ROOT::EnableThreadSafety();
#endif

	// set up numpy for usage in this module
	import_array();

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
	rpwa::py::exportAmpIntegralMatrixMetadata();
	rpwa::py::exportPhaseSpaceIntegral();
	rpwa::py::exportNBodyPhaseSpaceKinematics();
	rpwa::py::exportNBodyPhaseSpaceGenerator();
	rpwa::py::exportRandomNumberGenerator();
	rpwa::py::exportGeneratorManager();
	rpwa::py::exportGenerator();
	rpwa::py::exportGeneratorParameters();
	rpwa::py::exportGeneratorPickerFunctions();
	rpwa::py::exportBeamAndVertexGenerator();
#ifdef USE_BAT
	rpwa::py::exportImportanceSampler();
#endif
	rpwa::py::exportModelIntensity();
	rpwa::py::exportComplexMatrix();
	rpwa::py::exportFitResult();
	rpwa::py::exportPartialWaveFitHelper();
	rpwa::py::exportPhysUtils();
	rpwa::py::exportReportingUtilsEnvironment();
	rpwa::py::exportEventFileWriter();
	rpwa::py::exportEventMetadata();
	rpwa::py::exportHashCalculator();
	rpwa::py::exportAmplitudeFileWriter();
	rpwa::py::exportAmplitudeMetadata();
	rpwa::py::exportCalcAmplitude();
	rpwa::py::exportPwaLikelihood();
	rpwa::py::exportPwaFit();
	rpwa::py::exportGetMassShapes();
#ifdef USE_NLOPT
	rpwa::py::exportPwaNloptFit();
#endif

}
