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


set(PYUTILS_SUBDIR pyUtils)

set(DECAYAMPLITUDE_SUBDIR decayAmplitude)
set(GENERATORS_SUBDIR generators)
set(HIGHLEVELINTERFACE_SUBDIR highLevelInterface)
set(NBODYPHASESPACE_SUBDIR nBodyPhaseSpace)
set(PARTICLEDATA_SUBDIR particleData)
set(PARTIALWAVEFIT_SUBDIR partialWaveFit)
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
	${PYTHON_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	)


# source files that are compiled into library
set(SOURCES
	rootPwaPy.cc
	${PYUTILS_SUBDIR}/stlContainers_py.cc
	${PYUTILS_SUBDIR}/rootConverters_py.cc
	${DECAYAMPLITUDE_SUBDIR}/interactionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/fsVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/massDependence_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/productionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/diffractiveDissVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/decayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarCanonicalAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarHelicityAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/waveDescription_py.cc
	${DECAYAMPLITUDE_SUBDIR}/ampIntegralMatrix_py.cc
	${DECAYAMPLITUDE_SUBDIR}/phaseSpaceIntegral_py.cc
	${GENERATORS_SUBDIR}/generatorManager_py.cc
	${GENERATORS_SUBDIR}/generator_py.cc
	${GENERATORS_SUBDIR}/generatorParameters_py.cc
	${GENERATORS_SUBDIR}/generatorPickerFunctions_py.cc
	${GENERATORS_SUBDIR}/beamAndVertexGenerator_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/calcAmplitude_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/pwaFit_py.cc
	${NBODYPHASESPACE_SUBDIR}/nBodyPhaseSpaceGen_py.cc
	${NBODYPHASESPACE_SUBDIR}/randomNumberGenerator_py.cc
	${PARTIALWAVEFIT_SUBDIR}/complexMatrix_py.cc
	${PARTIALWAVEFIT_SUBDIR}/fitResult_py.cc
	${PARTICLEDATA_SUBDIR}/particleProperties_py.cc
	${PARTICLEDATA_SUBDIR}/particleDataTable_py.cc
	${PARTICLEDATA_SUBDIR}/particle_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeMetadata_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeTreeLeaf_py.cc
	${STORAGEFORMATS_SUBDIR}/eventFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/eventMetadata_py.cc
	${UTILITIES_SUBDIR}/utils_py.cc
	)
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


# link the shared library with the python bindings
add_custom_target(
	libRootPwaPyLink ALL
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}RootPwaPy${CMAKE_SHARED_LIBRARY_SUFFIX}	${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/libRootPwaPy.so
)


# produce python byte-code and move it to the build directory
add_custom_target(
	pyRootPwaPackageLinks ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/pyInterface/python_compiler.py ${CMAKE_SOURCE_DIR}/pyInterface/package ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa
)


add_dependencies(libRootPwaPyLink RootPwaPy pyRootPwaPackageLinks)


add_test(testPyRootPwa ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/pyInterface/testPyRootPwa.py)
