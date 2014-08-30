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
#//      build file for utilities
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# set include directories
set(INCLUDE_DIRECTORIES
	${UTILITIES_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	dataMetadata.cc
	)

# ROOT dictionary
set(ROOTPWAMETA_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaUtilsDict.cc)
set(ROOTPWAMETA_DICT_INCLUDE_DIRS
	${CMAKE_CURRENT_SOURCE_DIR}
)
root_generate_dictionary(
	"${ROOTPWAMETA_DICTIONARY}"
	"${ROOTPWAMETA_DICT_INCLUDE_DIRS}"
	"dataMetadata.h"
	"rootPwaMetaLinkDef.h"
)
set(SOURCES ${SOURCES} ${ROOTPWAMETA_DICTIONARY})


# library
set(THIS_LIB "RootPwaMeta")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	)
add_dependencies(${THIS_LIB} environment)

make_executable(convertDataFile convertDataFile.cc ${THIS_LIB} ${UTILITIES_LIBS})
make_executable(printMetadata   printMetadata.cc   ${THIS_LIB} ${UTILITIES_LIBS})
