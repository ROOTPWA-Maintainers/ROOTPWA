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


# set include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_NBODYPHASESPACE_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${PYTHON_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	)


# source files that are compiled into library
set(SOURCES
	interactionVertex.cc
	fsVertex.cc
	decayTopology.cc
	productionVertex.cc
	diffractiveDissVertex.cc
	leptoProductionVertex.cc
	nonInteractionVertex.cc
	isobarDecayVertex.cc
	isobarDecayTopology.cc
	massDependence.cc
	waveDescription.cc
	isobarAmplitude.cc
	isobarHelicityAmplitude.cc
	isobarCanonicalAmplitude.cc
	amplitudeTreeLeaf.cc
	amplitudeTreeHelper.cc
	ampIntegralMatrix.cc
	waveSetGenerator.cc
	phaseSpaceIntegral.cc
	)


# library name
lib_name_for_current_dir(THIS_LIB)


# ROOT dictionary
set(ROOTPWAAMP_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaAmpDict.cc)
root_generate_dictionary(
	${ROOTPWAAMP_DICTIONARY}
	waveDescription.h amplitudeTreeLeaf.h ampIntegralMatrix.h
	MODULE ${THIS_LIB}
	LINKDEF rootPwaAmpLinkDef.h
	)
set(SOURCES ${SOURCES} ${ROOTPWAAMP_DICTIONARY})


# library
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_NBODYPHASESPACE_LIB}"
	)


# executables
make_executable(calcIntegrals      calcIntegrals.cc      ${THIS_LIB} "${RPWA_UTILITIES_LIB}")
make_executable(ampDiff            ampDiff.cc            ${THIS_LIB} "${RPWA_UTILITIES_LIB}")
make_executable(createGraphDiagram createGraphDiagram.cc ${THIS_LIB} "${RPWA_UTILITIES_LIB}")
make_executable(generateWaveSet    generateWaveSet.cc    ${THIS_LIB} "${RPWA_UTILITIES_LIB}")


# build tests
add_subdirectory(test)
