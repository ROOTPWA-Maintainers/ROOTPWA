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
#//      build file for pwa library and amplitude executables
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'gamp' directory.")


set(CMAKE_CXX_FLAGS   "-O2 -Wall -Werror")
set(CMAKE_CXX_LDFLAGS "-O2")


# set include directories
set(INCLUDE_DIRECTORIES
	${UTILITIES_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/pwa2000/libpp
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	integral.cc
	keyParse.cc
	keyfile.cc
	keyScan.cc
	wave.cc
	Tgamp.cc
	)


# library
set(THIS_LIB "pwa")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	pp
	)


# executables
make_executable(gamp gamp.cc ${THIS_LIB})
make_executable(vamp vamp.cc ${THIS_LIB})
make_executable(int  int.cc  ${THIS_LIB})
