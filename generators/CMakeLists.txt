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
	SYSTEM ${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${TOOLS_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})

set(SOURCES_NBPSG
	nBodyPhaseSpaceGen.cc
	randomNumberGenerator.cc
	)

# ROOT dictionary
set(ROOTPWAGEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaGenDict.cc)
set(ROOTPWAGEN_DICT_INCLUDE_DIRS
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
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
set(SOURCES_NBPSG ${SOURCES_NBPSG} ${ROOTPWAGEN_DICTIONARY})

# library
set(NBPSG_LIB "RootPwaNBPSGen")
make_shared_library(
	"${NBPSG_LIB}"
	"${SOURCES_NBPSG}"
	"${ROOT_LIBS}"
	)


# source files that are compiled into library
set(SOURCES
	generator.cc
	generatorManager.cc
	diffractivePhaseSpace.cc
	beamAndVertexGenerator.cc
	generatorPickerFunctions.cc
	)

# library
set(THIS_LIB "RootPwaGen")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${ROOTPWA_LIBS}"
	"${ROOTPWATOOLS_LIBS}"
	"${NBPSG_LIB}"
	)


# executables
make_executable(generatePhaseSpace       generatePhaseSpace.cc           ${THIS_LIB} "${Libconfig_LIBS}" "${UTILITIES_LIBS}")
make_executable(weightEvents             weightEvents.cc                 ${THIS_LIB} "${AMPLITUDE_LIBS}" "${UTILITIES_LIBS}")

# build tests
add_subdirectory(test)
