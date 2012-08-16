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
#//      build file for RootPwaGen
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'generators' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${PWA2000_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${TOOLS_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	nBodyPhaseSpaceGen.cc
	TProductionAmp.cc
	TBWProductionAmp.cc
	TPWWeight.cc
	TDiffractivePhaseSpace.cc
	TPrimaryVertexGen.cc
	plotWeightedEvts.C  # this helper is included to the library to make it available to the gui
	)


# ROOT dictionary
set(ROOTPWAGEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaGenDict.cc)
set(ROOTPWAGEN_DICT_INCLUDE_DIRS
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${PWA2000_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${TOOLS_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
root_generate_dictionary(
	"${ROOTPWAGEN_DICTIONARY}"
	"${ROOTPWAGEN_DICT_INCLUDE_DIRS}"
	"nBodyPhaseSpaceGen.h"
	"rootPwaGenLinkDef.h"
	)
set(SOURCES ${SOURCES} ${ROOTPWAGEN_DICTIONARY})


# library
set(THIS_LIB "RootPwaGen")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${PWA2000_LIBS}"
	"${ROOTPWA_LIBS}"
	"${ROOTPWATOOLS_LIBS}"
	)


# executables
make_executable(genpw       genpw.cc           ${THIS_LIB} "${Libconfig_LIBS}")
make_executable(evtweight   evtweight.cc       ${THIS_LIB} "${AMPLITUDE_LIBS}")
make_executable(plotBWAmp   plotBWAmp.cc       ${THIS_LIB})
make_executable(doPlotWEvts plotWeightedEvts.C ${THIS_LIB})
