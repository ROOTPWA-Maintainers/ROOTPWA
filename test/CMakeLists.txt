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
#//      build file for tests
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()

add_subdirectory(decayAmplitude)
add_subdirectory(generators)
add_subdirectory(partialWaveFit)
add_subdirectory(utilities)


# set include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${RPWA_PARTIALWAVEFIT_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${RPWA_STORAGEFORMATS_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	)


# executables
make_executable(compareFitResultsIntegrals compareFitResultsIntegrals.cc "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_PARTIALWAVEFIT_LIB}" "${RPWA_UTILITIES_LIB}")


# test functionality of 'pyRootPwa' module
add_test(
	NAME testPyRootPwa
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testPyRootPwa.py
)
set_property(
	TEST testPyRootPwa
	APPEND
	PROPERTY ENVIRONMENT "PYTHONPATH=${PYTHON_LIBRARY_OUTPUT_PATH}:$ENV{PYTHONPATH}" "ROOTPWA=${CMAKE_SOURCE_DIR}"
)
