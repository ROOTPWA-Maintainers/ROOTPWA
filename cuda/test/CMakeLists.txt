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
#//      build file for CUDA tests
#//
#//
#// Author List:
#//      Boris Grube    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# set include directories
set(INCLUDE_DIRECTORIES
	SYSTEM ${Boost_INCLUDE_DIRS}
	SYSTEM ${CUDA_SAMPLES_INC_DIR}
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${CUDACODE_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# compile HOST executables
make_executable(testLikelihoodMockup testLikelihoodMockup.cc ${CUDA_CUTIL_LIBRARY} "${CUDACODE_LIBS}" "${ROOT_LIBS}")

cuda_add_executable(testComplexMemoryBandwidth testComplexMemoryBandwidth.cu)
target_link_libraries(testComplexMemoryBandwidth ${CUDA_CUTIL_LIBRARY})

cuda_add_executable(benchmarkMemoryBandwidth benchmarkMemoryBandwidth.cu)
target_link_libraries(benchmarkMemoryBandwidth ${CUDA_CUTIL_LIBRARY})
