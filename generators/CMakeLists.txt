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
#//      build file for generators library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# set include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${RPWA_NBODYPHASESPACE_INCLUDE_DIR}
	${RPWA_PARTIALWAVEFIT_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	${RPWA_STORAGEFORMATS_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
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
lib_name_for_current_dir(THIS_LIB)
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_PARTIALWAVEFIT_LIB}"
	"${RPWA_NBODYPHASESPACE_LIB}"
	"${RPWA_STORAGEFORMATS_LIB}"
	)


# executables
make_executable(generatePhaseSpace generatePhaseSpace.cc ${THIS_LIB} "${Libconfig_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(weightEvents       weightEvents.cc       ${THIS_LIB} "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_UTILITIES_LIB}")


# build tests
add_subdirectory(test)
