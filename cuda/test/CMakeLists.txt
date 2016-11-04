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


# adds CUDA executable
# additional libraries that should be linked to can be given as optional arguments
function(make_cuda_executable EXE_NAME SOURCES)
        message(STATUS ">>> setting up CUDA executable '${EXE_NAME}'")
        if(DEBUG_OUTPUT)
                message(STATUS "make_cuda_executable was called with the following arguments:
        EXE_NAME = '${EXE_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
        endif()
        cuda_add_executable(${EXE_NAME} ${SOURCES})
        # proccess link libraries in additional arguments
        foreach(_LIB ${ARGN})
                target_link_libraries(${EXE_NAME} ${_LIB})
        endforeach()
        unset(_LIB)
endfunction(make_cuda_executable)


# set include directories
include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_CUDA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${CUDA_SAMPLES_INC_DIR}
	${ROOT_INCLUDE_DIR}
	)


# compile HOST executables
make_executable(testLikelihoodMockup testLikelihoodMockup.cc ${CUDA_CUTIL_LIBRARY} "${RPWA_CUDA_LIB}" "${ROOT_LIBS}")


# compile CUDA executables
make_cuda_executable(testComplexMemoryBandwidth testComplexMemoryBandwidth.cu ${CUDA_CUTIL_LIBRARY})
make_cuda_executable(benchmarkMemoryBandwidth   benchmarkMemoryBandwidth.cu   ${CUDA_CUTIL_LIBRARY})
