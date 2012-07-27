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
#// File and Version Information:
#// $Rev::                             $: revision of last commit
#// $Author::                          $: author of last commit
#// $Date::                            $: date of last commit
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


add_definitions(-Wno-deprecated)  # suppress deprecated header warnings from boost libraries when using gcc 4.4


# set include directories
set(INCLUDE_DIRECTORIES
	${PYTHON_INCLUDE_DIRS}
	${uprootedBoost_INCLUDE_DIRS}  # override some Boost libraries with local copies
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${PWA2000_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	particleProperties.cc
	particleDataTable.cc
	particle.cc
	interactionVertex.cc
	fsVertex.cc
	decayTopology.cc
	productionVertex.cc
	diffractiveDissVertex.cc
	leptoProductionVertex.cc
	isobarDecayVertex.cc
	isobarDecayTopology.cc
	massDependence.cc
	waveDescription.cc
	evtTreeHelper.cc
	isobarAmplitude.cc
	isobarHelicityAmplitude.cc
	isobarCanonicalAmplitude.cc
	amplitudeTreeLeaf.cc
	amplitudeTreeHelper.cc
	ampIntegralMatrix.cc
	waveSetGenerator.cc
	)

set(PYSOURCES
	python_stuff.cc
	)


# ROOT dictionary
set(ROOTPWAAMP_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaAmpDict.cc)
set(ROOTPWAAMP_DICT_INCLUDE_DIRS
	${CMAKE_CURRENT_SOURCE_DIR}
	)
root_generate_dictionary(
	"${ROOTPWAAMP_DICTIONARY}"
	"${ROOTPWAAMP_DICT_INCLUDE_DIRS}"
	"waveDescription.h;amplitudeTreeLeaf.h;ampIntegralMatrix.h"
	"rootPwaAmpLinkDef.h"
	)
set(SOURCES ${SOURCES} ${ROOTPWAAMP_DICTIONARY})


# library
set(THIS_LIB "RootPwaAmp")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${BOOST_LIBS}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${PWA2000_LIBS}"
	)

make_shared_library(
	"pylib"
	"${PYSOURCES}"
	"${PYTHON_LIBRARIES}"
	"${Boost_Libraries}"
	"${BOOST_LIBS}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${PWA2000_LIBS}"
	)


# executables
make_executable(calcAmplitudes     calcAmplitudes.cc     ${THIS_LIB})
make_executable(calcIntegrals      calcIntegrals.cc      ${THIS_LIB})
make_executable(checkKeyFile       checkKeyFile.cc       ${THIS_LIB})
make_executable(ampDiff            ampDiff.cc            ${THIS_LIB})
make_executable(createGraphDiagram createGraphDiagram.cc ${THIS_LIB})
make_executable(generateWaveSet    generateWaveSet.cc    ${THIS_LIB})

make_executable(testParticle         testParticle.cc         ${THIS_LIB})
make_executable(testDecayTopology    testDecayTopology.cc    ${THIS_LIB})
make_executable(testWaveDescription  testWaveDescription.cc  ${THIS_LIB})
make_executable(testWaveSetGenerator testWaveSetGenerator.cc ${THIS_LIB})
make_executable(testAmplitude        testAmplitude.cc        ${THIS_LIB})
make_executable(testAmplitudeTree    testAmplitudeTree.cc    ${THIS_LIB})
make_executable(testIntegral         testIntegral.cc         ${THIS_LIB})
