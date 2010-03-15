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
#//      main build file for ROOTPWA
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)


# set verbosity
SET(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
                               # the same effect can be achieved by calling 'make VERBOSE=1'


# define project
PROJECT(rootpwa)
ENABLE_TESTING()


# define build types
# set a default build type for single-configuration CMake generators, if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_FLAGS "-Wall -Werror")  # common flags
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_LDFLAGS_DEBUG "-g")


# define output directories
SET(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


# setup Boost
# environment variable $BOOST_ROOT is expected to point to non-standard locations
SET(Boost_INCLUDE_DIRS $ENV{BOOST_ROOT})
IF(Boost_INCLUDE_DIRS)  # ugly CMAKE 2.4 workaround
  MESSAGE(STATUS "Using Boost include directory '${Boost_INCLUDE_DIRS}'")
ELSE(Boost_INCLUDE_DIRS)
  FIND_PACKAGE(Boost)
  IF(Boost_FOUND)
    MESSAGE(STATUS "Using Boost include directory '${Boost_INCLUDE_DIRS}'")
    MESSAGE(STATUS "Using Boost library directory '${Boost_LIBRARY_DIRS}'")
    FIND_LIBRARY(BOOST_GRAPH_LIB
      libboost_graph.so
      PATHS ${Boost_LIBRARY_DIRS}
    )
    IF(${BOOST_GRAPH_LIB} MATCHES "BOOST_GRAPH_LIB-NOTFOUND")
      MESSAGE(STATUS "Could not find 'libboost_graph.so'. Disabling some features.")
      SET(BOOST_GRAPH_LIB "")
    ELSE(${BOOST_GRAPH_LIB} MATCHES "BOOST_GRAPH_LIB-NOTFOUND")
      MESSAGE(STATUS "Using Boost Graph Library '${BOOST_GRAPH_LIB}'")
    ENDIF(${BOOST_GRAPH_LIB} MATCHES "BOOST_GRAPH_LIB-NOTFOUND")
  ELSE(Boost_FOUND)
    MESSAGE(FATAL_ERROR "Could not find Boost installation. Is environment variable BOOST_ROOT=${BOOST_ROOT} set correctly?")
  ENDIF(Boost_FOUND)
ENDIF(Boost_INCLUDE_DIRS)
IF(BOOST_GRAPH_LIB)
  ADD_DEFINITIONS(-D'BOOST_GRAPH_LIB=1')
ELSE(BOOST_GRAPH_LIB)
  ADD_DEFINITIONS(-D'BOOST_GRAPH_LIB=0')
ENDIF(BOOST_GRAPH_LIB)
SET(BOOST_LIBS
  ${BOOST_GRAPH_LIB}
)


# make SVN version string accessible via predefined macro SVN_VERSION
FIND_PACKAGE(Subversion)
IF(Subversion_FOUND)
  # unfortunately CMAKE only parses 'svn info'
  FIND_PROGRAM(SVNVERSION_EXECUTABLE
    svnversion
  )
  IF(${SVNVERSION_EXECUTABLE} MATCHES "SVNVERSION_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "Could not find subversion command 'svnversion'. Repository version unknown.")
  ELSE(${SVNVERSION_EXECUTABLE} MATCHES "SVNVERSION_EXECUTABLE-NOTFOUND")
    EXEC_PROGRAM(${SVNVERSION_EXECUTABLE}
      ARGS "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE SVN_VERSION
      RETURN_VALUE SVNVERSION_RETURN
    )
    IF(NOT ${SVNVERSION_RETURN})
      MESSAGE(STATUS "Subversion repository version is '${SVN_VERSION}'")
    ELSE(NOT ${SVNVERSION_RETURN})
      MESSAGE(STATUS "Error running 'svnversion'. Repository version unknown.")
      SET(SVN_VERSION "")
    ENDIF(NOT ${SVNVERSION_RETURN})
  ENDIF(${SVNVERSION_EXECUTABLE} MATCHES "SVNVERSION_EXECUTABLE-NOTFOUND")
ELSE(Subversion_FOUND)
  MESSAGE(STATUS "Could not find subversion installation. Repository version unknown.")
ENDIF(Subversion_FOUND)
IF(SVN_VERSION)
  ADD_DEFINITIONS(-D'SVN_VERSION=\"${SVN_VERSION}\"')
ELSE(SVN_VERSION)
  ADD_DEFINITIONS(-D'SVN_VERSION=\"\"')
ENDIF(SVN_VERSION)


# setup ROOT includes and libraries
SET(ROOT_CONFIG_SEARCHPATH $ENV{ROOTSYS}/bin)
FIND_PROGRAM(ROOT_CONFIG_EXECUTABLE
  root-config
  PATHS ${ROOT_CONFIG_SEARCHPATH}
  NO_DEFAULT_PATH
)
IF(${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
  MESSAGE(FATAL_ERROR "ROOT not installed in the searchpath and ROOTSYS is not set. Please set ROOTSYS or add the path to your ROOT installation in the Macro FindROOT.cmake in th subdirectory cmake/modules.")
ELSE(${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
  STRING(REGEX REPLACE "(^.*)/bin/root-config" "\\1" test ${ROOT_CONFIG_EXECUTABLE})
  SET($ENV{ROOTSYS} ${test})
ENDIF(${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND") 
# ROOT version
EXEC_PROGRAM(${ROOT_CONFIG_EXECUTABLE}
  ARGS "--version --svn-revision"
  OUTPUT_VARIABLE ROOT_VERSION
)
MESSAGE(STATUS "Linking against ROOT version ${ROOT_VERSION}")
# ROOT includes
EXEC_PROGRAM(${ROOT_CONFIG_EXECUTABLE}
  ARGS "--incdir"
  OUTPUT_VARIABLE root_headers
)
SET(ROOT_INCLUDE_DIR ${root_headers})
MESSAGE(STATUS "Using ROOT include directory '${ROOT_INCLUDE_DIR}'")
# ROOT libraries
EXEC_PROGRAM(${ROOT_CONFIG_EXECUTABLE}
  ARGS "--libdir"
  OUTPUT_VARIABLE root_libs
)
SET(ROOT_LIB ${root_libs})
MESSAGE(STATUS "Linking against ROOT libraries in '${ROOT_LIB}'")


# Recurse into the given subdirectories. This does not actually cause
# another cmake executable to run. The same process will walk through
# the project's entire directory structure.
ADD_SUBDIRECTORY(pwa2000)
ADD_SUBDIRECTORY(tools)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(generators)
ADD_SUBDIRECTORY(keygen)
ADD_SUBDIRECTORY(amplitude)
#ADD_SUBDIRECTORY(gui)
 
