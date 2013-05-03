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
#//      build file for RootPwaTools library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'regge' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${PWA2000_INCLUDE_DIR}
	${TOOLS_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	reggeprop.cc
	pipiamp.cc
	pipiphaseshift.cc
	)


# ROOT dictionary
#set(REGGE_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaReggeDict.cc)
#root_generate_dictionary(
#  "${TOOLS_DICTIONARY}"
##  "${CMAKE_CURRENT_SOURCE_DIR}"
#  "rootPwaReggeLinkDef.h"
#)
#set(SOURCES ${SOURCES} ${TOOLS_DICTIONARY})


# library
set(THIS_LIB "RootPwaRegge")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${PWA2000_LIBS}"
	"${ROOTPWATOOLS_LIBS}"
	)


# executables
make_executable(testRegge testRegge.cc ${THIS_LIB})
