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
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${RPWA_GENERATORS_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	)


# executables
make_executable(testParticle                testParticle.cc                  "${RPWA_DECAYAMPLITUDE_LIB}")
make_executable(testDecayTopology           testDecayTopology.cc             "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_UTILITIES_LIB}")
make_executable(testWaveDescription         testWaveDescription.cc           "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_UTILITIES_LIB}")
make_executable(testWaveSetGenerator        testWaveSetGenerator.cc          "${RPWA_DECAYAMPLITUDE_LIB}")
make_executable(testAmplitudeTree           testAmplitudeTree.cc             "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_UTILITIES_LIB}")
make_executable(testIntegral                testIntegral.cc                  "${RPWA_DECAYAMPLITUDE_LIB}")
make_executable(testPeakShapeInstableIsobar testPeakShapeInstableIsobar.cc   "${RPWA_DECAYAMPLITUDE_LIB}" "${RPWA_GENERATORS_LIB}")
make_executable(testPeakShapeInstableIsobar2 testPeakShapeInstableIsobar2.cc "${RPWA_DECAYAMPLITUDE_LIB}")
