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
#//      build file for RootPwa library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'nlopt' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${NLopt_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})

# executables
make_executable(pwaNloptFit   pwaNloptFit.cc       ${ROOTPWA_LIBS} ${NLopt_LIBS} ${UTILITIES_LIBS})
