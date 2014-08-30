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


set(AMPLITUDE_SUBDIR amplitude)
set(GENERATORS_SUBDIR generators)
set(METADATA_SUBDIR metadata)
set(PYUTILS_SUBDIR pyUtils)
set(SRC_SUBDIR src)
set(UTILITIES_SUBDIR utilities)



# set include directories
set(INCLUDE_DIRECTORIES
	${PYTHON_INCLUDE_DIRS}
	SYSTEM ${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${GEN_INCLUDE_DIR}
	${METADATA_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${AMPLITUDE_SUBDIR}
	${GENERATORS_SUBDIR}
	${METADATA_SUBDIR}
	${PYUTILS_SUBDIR}
	${SRC_SUBDIR}
	${UTILITIES_SUBDIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	rootPwaPy.cc
	${PYUTILS_SUBDIR}/stlContainers_py.cc
	${PYUTILS_SUBDIR}/rootConverters_py.cc
	${AMPLITUDE_SUBDIR}/particleProperties_py.cc
	${AMPLITUDE_SUBDIR}/particleDataTable_py.cc
	${AMPLITUDE_SUBDIR}/particle_py.cc
	${AMPLITUDE_SUBDIR}/interactionVertex_py.cc
	${AMPLITUDE_SUBDIR}/fsVertex_py.cc
	${AMPLITUDE_SUBDIR}/massDependence_py.cc
	${AMPLITUDE_SUBDIR}/isobarDecayVertex_py.cc
	${AMPLITUDE_SUBDIR}/productionVertex_py.cc
	${AMPLITUDE_SUBDIR}/diffractiveDissVertex_py.cc
	${AMPLITUDE_SUBDIR}/decayTopology_py.cc
	${AMPLITUDE_SUBDIR}/isobarDecayTopology_py.cc
	${AMPLITUDE_SUBDIR}/isobarAmplitude_py.cc
	${AMPLITUDE_SUBDIR}/isobarCanonicalAmplitude_py.cc
	${AMPLITUDE_SUBDIR}/isobarHelicityAmplitude_py.cc
	${AMPLITUDE_SUBDIR}/waveDescription_py.cc
	${AMPLITUDE_SUBDIR}/amplitudeTreeLeaf_py.cc
	${AMPLITUDE_SUBDIR}/ampIntegralMatrix_py.cc
	${AMPLITUDE_SUBDIR}/phaseSpaceIntegral_py.cc
	${GENERATORS_SUBDIR}/nBodyPhaseSpaceGen_py.cc
	${GENERATORS_SUBDIR}/randomNumberGenerator_py.cc
	${GENERATORS_SUBDIR}/generatorManager_py.cc
	${GENERATORS_SUBDIR}/generator_py.cc
	${GENERATORS_SUBDIR}/generatorParameters_py.cc
	${GENERATORS_SUBDIR}/generatorPickerFunctions_py.cc
	${GENERATORS_SUBDIR}/beamAndVertexGenerator_py.cc
	${METADATA_SUBDIR}/dataFileWriter_py.cc
	${SRC_SUBDIR}/complexMatrix_py.cc
	${SRC_SUBDIR}/fitResult_py.cc
	${UTILITIES_SUBDIR}/utils_py.cc
	)


# library
set(THIS_LIB "RootPwaPy")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${PYTHON_LIBRARIES}"
	"${Boost_PYTHON_LIBRARY}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${AMPLITUDE_LIBS}"
	"${ROOTPWAGEN_LIBS}"
	"${UTILITIES_LIBS}"
	"${METADATA_LIBS}"
	)

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
