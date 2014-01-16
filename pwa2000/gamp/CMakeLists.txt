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


# override compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
# make sure flags are supported by compiler
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-unneeded-internal-declaration COMPILER_HAS_FLAG_NO_UNNEEDED_INTERNAL_DECLARATION)
if(COMPILER_HAS_FLAG_NO_UNNEEDED_INTERNAL_DECLARATION)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unneeded-internal-declaration")
endif()
unset(COMPILER_HAS_FLAG_NO_UNNEEDED_INTERNAL_DECLARATION)
check_cxx_compiler_flag(-Wno-unused-variable COMPILER_HAS_FLAG_NO_UNUSED_VARIABLE)
if(COMPILER_HAS_FLAG_NO_UNUSED_VARIABLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
endif()
unset(COMPILER_HAS_FLAG_NO_UNUSED_VARIABLE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")


# set include directories
set(INCLUDE_DIRECTORIES
	${UTILITIES_INCLUDE_DIR}
	SYSTEM ${Boost_INCLUDE_DIRS}
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
