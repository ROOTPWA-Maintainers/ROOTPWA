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


message(STATUS "")
message(STATUS ">>> Setting up 'pyInterface' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${PYTHON_INCLUDE_DIRS}
	SYSTEM ${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
#	${PWA2000_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	rootPwaPy.cc
	stlContainers_py.cc
	rootConverters_py.cc
	particleProperties_py.cc
	particleDataTable_py.cc
	particle_py.cc
	interactionVertex_py.cc
	fsVertex_py.cc
	massDependence_py.cc
	isobarDecayVertex_py.cc
	productionVertex_py.cc
	diffractiveDissVertex_py.cc
	decayTopology_py.cc
	isobarDecayTopology_py.cc
	isobarAmplitude_py.cc
	isobarCanonicalAmplitude_py.cc
	isobarHelicityAmplitude_py.cc
	waveDescription_py.cc
	amplitudeTreeLeaf_py.cc
	ampIntegralMatrix_py.cc
	utils_py.cc
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
	"${PWA2000_LIBS}"
	"${AMPLITUDE_LIBS}"
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
