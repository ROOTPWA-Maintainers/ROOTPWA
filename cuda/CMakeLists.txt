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
#//      test build file for CUDA
#//
#//
#// Author List:
#//      Boris Grube    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


if(USE_CUDA)

  set(CMAKE_CXX_FLAGS "-O3")
  message(STATUS "Using CXX flags '${CMAKE_CXX_FLAGS}' for CUDA compilation")

  # define for which shader models
  #set(_SHADER_MODELS "13" "20")
  set(_SHADER_MODELS "13")
  foreach(_SHADER_MODEL ${_SHADER_MODELS})
    list(APPEND CUDA_NVCC_FLAGS "--generate-code=arch=compute_${_SHADER_MODEL},code=\\\"sm_${_SHADER_MODEL},compute_${_SHADER_MODEL}\\\"")
  endforeach()
  list(APPEND CUDA_NVCC_FLAGS "--ptxas-options=-v")
  message(STATUS "Compiling using CUDA nvcc flags '${CUDA_NVCC_FLAGS}'")

  # set include directories
  set(INCLUDE_DIRECTORIES
    ${Boost_INCLUDE_DIRS}
    ${CUDA_SDK_ROOT_DIR}/common/inc
		${ROOT_INCLUDE_DIR}
    ${UTILITIES_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  include_directories(${INCLUDE_DIRECTORIES})

  # find cuda utility library
  find_library(CUDA_CUTIL_LIBRARY
    NAMES cutil_${CMAKE_SYSTEM_PROCESSOR}
    PATHS ${CUDA_SDK_ROOT_DIR}/lib
    NO_DEFAULT_PATH
    )
  if(NOT CUDA_CUTIL_LIBRARY)
    message(FATAL_ERROR "Could not find CUDA utility library 'cutil_${CMAKE_SYSTEM_PROCESSOR}' in '${CUDA_SDK_ROOT_DIR}/lib'")
  else()
    message(STATUS "Using CUDA utility library '${CUDA_CUTIL_LIBRARY}'")
  endif()

  # compile CUDA library
  cuda_add_library(cuda likelihoodInterface.cu SHARED)

  # compile CUDA executable
  make_executable(testLikelihoodMockup testLikelihoodMockup.cc)
  target_link_libraries(testLikelihoodMockup ${CUDA_CUTIL_LIBRARY} cuda ${ROOT_LIBS})

  cuda_add_executable(testComplexMemoryBandwidth testComplexMemoryBandwidth.cu)
  target_link_libraries(testComplexMemoryBandwidth ${CUDA_CUTIL_LIBRARY})

  cuda_add_executable(benchmarkMemoryBandwidth benchmarkMemoryBandwidth.cu)
  target_link_libraries(benchmarkMemoryBandwidth ${CUDA_CUTIL_LIBRARY})

endif()
