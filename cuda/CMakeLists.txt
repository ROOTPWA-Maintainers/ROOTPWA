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


message(STATUS "")
message(STATUS ">>> Setting up 'cuda' directory.")


# -Wall reports warnings in system libraries; so we cannot use -Werror
set(CMAKE_CXX_FLAGS "-Wall")
message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS}' for CUDA host code")
message(STATUS "Using CUDA nvcc flags '${CUDA_NVCC_FLAGS}'")


# set include directories
set(INCLUDE_DIRECTORIES
	SYSTEM ${Boost_INCLUDE_DIRS}
	SYSTEM ${CUDA_SAMPLES_INC_DIR}
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
  ${AMPLITUDE_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# compile CUDA library
cuda_add_library(RootPwaCuda    likelihoodInterface.cu  SHARED)
cuda_add_library(RootPwaCudaAmp helampltestInterface.cu SHARED)


# compile HOST executables
make_executable(testLikelihoodMockup testLikelihoodMockup.cc)
target_link_libraries(testLikelihoodMockup RootPwaCuda ${ROOT_LIBS})

make_executable(helampltest helampltest.cc)
target_link_libraries(helampltest ${ROOT_LIBS} ${AMPLITUDE_LIBS})

make_executable(testAmplitudeSpeed testAmplitudeSpeed.cc)
target_link_libraries(testAmplitudeSpeed ${ROOT_LIBS} ${AMPLITUDE_LIBS} RootPwaCudaAmp RootPwaCuda)


# compile DEVICE executables
cuda_add_executable(testComplexMemoryBandwidth testComplexMemoryBandwidth.cu)
#target_link_libraries(testComplexMemoryBandwidth)

cuda_add_executable(benchmarkMemoryBandwidth benchmarkMemoryBandwidth.cu)
#target_link_libraries(benchmarkMemoryBandwidth)

cuda_add_executable(square square.cu)
target_link_libraries(square ${ROOT_LIBS})

cuda_add_executable(templtest templtest.cu)
target_link_libraries(templtest ${ROOT_LIBS})

cuda_add_executable(bwig bwig.cu)
target_link_libraries(bwig ${ROOT_LIBS})

cuda_add_executable(treesum treesum.cu)
target_link_libraries(treesum ${ROOT_LIBS})
