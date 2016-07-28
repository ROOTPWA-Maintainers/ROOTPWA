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
#//      build file for nBodyPhaseSpace library
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
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	)


# source files that are compiled into library
set(SOURCES
	nBodyPhaseSpaceGenerator.cc
	nBodyPhaseSpaceKinematics.cc
	randomNumberGenerator.cc
	)


# library name
lib_name_for_current_dir(THIS_LIB)


# ROOT dictionary
set(ROOTPWAGEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/dict.cc)
root_generate_dictionary(
	${ROOTPWAGEN_DICTIONARY}
	nBodyPhaseSpaceGenerator.h nBodyPhaseSpaceKinematics.h
	MODULE ${THIS_LIB}
	LINKDEF linkdef.h
	)
set(SOURCES ${SOURCES} ${ROOTPWAGEN_DICTIONARY})


# library
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	)
