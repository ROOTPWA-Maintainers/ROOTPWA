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
#//      build file for amplitude library tests
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${GEN_INCLUDE_DIR}
	${PARTICLEDATA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# executables
make_executable(testParticle                testParticle.cc                  "${AMPLITUDE_LIBS}" "${PARTICLEDATA_LIBS}")
make_executable(testDecayTopology           testDecayTopology.cc             "${AMPLITUDE_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(testWaveDescription         testWaveDescription.cc           "${AMPLITUDE_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(testWaveSetGenerator        testWaveSetGenerator.cc          "${AMPLITUDE_LIBS}")
make_executable(testAmplitude               testAmplitude.cc                 "${AMPLITUDE_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(testAmplitudeTree           testAmplitudeTree.cc             "${AMPLITUDE_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(testIsospinSym              testIsospinSym.cc                "${AMPLITUDE_LIBS}" "${RPWA_UTILITIES_LIB}")
make_executable(testIntegral                testIntegral.cc                  "${AMPLITUDE_LIBS}")
make_executable(testPeakShapeInstableIsobar testPeakShapeInstableIsobar.cc   "${AMPLITUDE_LIBS}" "${ROOTPWAGEN_LIBS}")
make_executable(testPeakShapeInstableIsobar2 testPeakShapeInstableIsobar2.cc "${AMPLITUDE_LIBS}")
