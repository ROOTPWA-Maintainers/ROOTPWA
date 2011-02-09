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
#// $Rev:: 496                         $: revision of last commit
#// $Author:: bgrube                   $: author of last commit
#// $Date:: 2010-11-19 11:04:11 +0100 #$: date of last commit
#//
#// Description:
#//      test build file for CUDA
#//
#//
#// Author List:
#//      Boris Grube    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "Building MPI code")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LDFLAGS "${CMAKE_CXX_LDFLAGS} ${MPI_LINK_FLAGS}")
message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS}' for MPI code")
message(STATUS "Using CXX linker flags '${CMAKE_CXX_LDFLAGS}' for MPI code")


# set include directories
set(INCLUDE_DIRECTORIES
	${MPI_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )
include_directories(${INCLUDE_DIRECTORIES})


# executables
make_executable(testMpiPerformance.cpp testMpiPerformance.cpp ${MPI_LIBRARIES} ${Boost_MPI_LIBS})
