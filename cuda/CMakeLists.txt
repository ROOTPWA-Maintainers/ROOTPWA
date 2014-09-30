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
#//      build file for CUDA code
#//
#//
#// Author List:
#//      Boris Grube    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# -Wall reports warnings in system libraries; so we cannot use -Werror
set(CMAKE_CXX_FLAGS "-Wall")
message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS}' for CUDA host code")
message(STATUS "Using CUDA nvcc flags '${CUDA_NVCC_FLAGS}'")


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${CUDA_SAMPLES_INC_DIR}
	${ROOT_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	likelihoodInterface.cu
	)


# compile CUDA library
set(_THIS_LIB "RootPwaCuda")
cuda_add_library(${_THIS_LIB} ${SOURCES} SHARED)


# build tests
add_subdirectory(test)
