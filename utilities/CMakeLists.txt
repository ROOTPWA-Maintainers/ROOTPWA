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
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	SYSTEM
	${ROOT_INCLUDE_DIR}
)


# produce new version of environment printing variables
ADD_CUSTOM_TARGET(
	environment
	${CMAKE_COMMAND} -D SRC=${CMAKE_CURRENT_SOURCE_DIR}/environment.h.in
	                 -D DST=${CMAKE_CURRENT_SOURCE_DIR}/environment.h
	                 -D CMAKE_HOST_SYSTEM_NAME=${CMAKE_HOST_SYSTEM_NAME}
	                 -D CMAKE_HOST_SYSTEM_PROCESSOR=${CMAKE_HOST_SYSTEM_PROCESSOR}
	                 -D CMAKE_HOST_SYSTEM_VERSION=${CMAKE_HOST_SYSTEM_VERSION}
	                 -D NMB_CPU_CORES=${NMB_CPU_CORES}
	                 -D HOSTNAME=${HOSTNAME}
	                 -D USER=${USER}
	                 -D CMAKE_SOURCE_DIR_=${CMAKE_SOURCE_DIR}
	                 -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	                 -D Boost_LIBRARY_VERSION=${Boost_LIBRARY_VERSION}
	                 -D Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}
	                 -D Libconfig_VERSION=${Libconfig_VERSION}
	                 -D Libconfig_DIR=${Libconfig_DIR}
	                 -D ROOTSYS=${ROOTSYS}
	                 -D CUDA_VERSION=${CUDA_VERSION}
	                 -D CUDA_LIB_DIRS=${CUDA_LIB_DIRS}
	                 -D PYTHONLIBS_VERSION_STRING=${PYTHONLIBS_VERSION_STRING}
	                 -D PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}
	                 -P ${CMAKE_CURRENT_SOURCE_DIR}/environment.cmake
)


# source files that are compiled into library
set(SOURCES
	reportingUtilsEnvironment.cc
	)


# library
lib_name_for_current_dir(THIS_LIB)
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	)
add_dependencies(${THIS_LIB} environment)


# build tests
add_subdirectory(test)
