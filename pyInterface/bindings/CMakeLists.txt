#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      build file for amplitude library
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


set(DECAYAMPLITUDE_SUBDIR decayAmplitude)
set(GENERATORS_SUBDIR generators)
set(HIGHLEVELINTERFACE_SUBDIR highLevelInterface)
set(NBODYPHASESPACE_SUBDIR nBodyPhaseSpace)
set(PARTICLEDATA_SUBDIR particleData)
set(PARTIALWAVEFIT_SUBDIR partialWaveFit)
set(PYUTILS_SUBDIR pyUtils)
set(STORAGEFORMATS_SUBDIR storageFormats)
set(UTILITIES_SUBDIR utilities)


# set include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${PYUTILS_SUBDIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${DECAYAMPLITUDE_SUBDIR}
	${RPWA_GENERATORS_INCLUDE_DIR}
	${GENERATORS_SUBDIR}
	${RPWA_HIGHLEVELINTERFACE_INCLUDE_DIR}
	${HIGHLEVELINTERFACE_SUBDIR}
	${RPWA_NBODYPHASESPACE_INCLUDE_DIR}
	${NBODYPHASESPACE_SUBDIR}
	${RPWA_PARTIALWAVEFIT_INCLUDE_DIR}
	${PARTIALWAVEFIT_SUBDIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${PARTICLEDATA_SUBDIR}
	${RPWA_STORAGEFORMATS_INCLUDE_DIR}
	${STORAGEFORMATS_SUBDIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	${UTILITIES_SUBDIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${NUMPY_INCLUDE_DIR}
	${PYTHON_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	)
if(USE_BAT)
	include_directories(
		SYSTEM
		${BAT_INCLUDE_DIR}
	)
endif()


# source files that are compiled into library
set(SOURCES
	rootPwaPy.cc
	${PYUTILS_SUBDIR}/rootConverters_py.cc
	${PYUTILS_SUBDIR}/stlContainers_py.cc
	${DECAYAMPLITUDE_SUBDIR}/ampIntegralMatrix_py.cc
	${DECAYAMPLITUDE_SUBDIR}/ampIntegralMatrixMetadata_py.cc
	${DECAYAMPLITUDE_SUBDIR}/decayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/diffractiveDissVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/fsVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/interactionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarCanonicalAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarHelicityAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/massDependence_py.cc
	${DECAYAMPLITUDE_SUBDIR}/phaseSpaceIntegral_py.cc
	${DECAYAMPLITUDE_SUBDIR}/productionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/waveDescription_py.cc
	${GENERATORS_SUBDIR}/beamAndVertexGenerator_py.cc
	${GENERATORS_SUBDIR}/generator_py.cc
	${GENERATORS_SUBDIR}/generatorManager_py.cc
	${GENERATORS_SUBDIR}/generatorParameters_py.cc
	${GENERATORS_SUBDIR}/generatorPickerFunctions_py.cc
	${GENERATORS_SUBDIR}/modelIntensity_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/calcAmplitude_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/pwaFit_py.cc
	${NBODYPHASESPACE_SUBDIR}/nBodyPhaseSpaceGenerator_py.cc
	${NBODYPHASESPACE_SUBDIR}/nBodyPhaseSpaceKinematics_py.cc
	${NBODYPHASESPACE_SUBDIR}/randomNumberGenerator_py.cc
	${PARTIALWAVEFIT_SUBDIR}/complexMatrix_py.cc
	${PARTIALWAVEFIT_SUBDIR}/fitResult_py.cc
	${PARTIALWAVEFIT_SUBDIR}/partialWaveFitHelper_py.cc
	${PARTIALWAVEFIT_SUBDIR}/pwaLikelihood_py.cc
	${PARTICLEDATA_SUBDIR}/particle_py.cc
	${PARTICLEDATA_SUBDIR}/particleDataTable_py.cc
	${PARTICLEDATA_SUBDIR}/particleProperties_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeMetadata_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeTreeLeaf_py.cc
	${STORAGEFORMATS_SUBDIR}/eventFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/eventMetadata_py.cc
	${STORAGEFORMATS_SUBDIR}/hashCalculator_py.cc
	${UTILITIES_SUBDIR}/physUtils_py.cc
	${UTILITIES_SUBDIR}/reportingUtilsEnvironment_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/getMassShapes_py.cc
	)
if(USE_BAT)
	LIST(APPEND SOURCES ${GENERATORS_SUBDIR}/importanceSampler_py.cc)
endif()
if(USE_NLOPT)
	LIST(APPEND SOURCES ${HIGHLEVELINTERFACE_SUBDIR}/pwaNloptFit_py.cc)
endif()


# library
set(THIS_LIB "RootPwaPy")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${PYTHON_LIBRARIES}"
	"${Boost_PYTHON_LIBRARY}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${RPWA_DECAYAMPLITUDE_LIB}"
	"${RPWA_GENERATORS_LIB}"
	"${RPWA_HIGHLEVELINTERFACE_LIB}"
	"${RPWA_NBODYPHASESPACE_LIB}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_STORAGEFORMATS_LIB}"
	"${RPWA_UTILITIES_LIB}"
	)
if(USE_NLOPT)
	target_link_libraries(${THIS_LIB} "${NLopt_LIBS}")
endif()
