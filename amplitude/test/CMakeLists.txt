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
#// $Rev:: 812                         $: revision of last commit
#// $Author:: bgrube                   $: author of last commit
#// $Date:: 2011-11-11 19:55:55 +0100 #$: date of last commit
#//
#// Description:
#//      build file for amplitude library 
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# set include directories
set(INCLUDE_DIRECTORIES
  ${Boost_INCLUDE_DIRS}
  ${Libconfig_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
  ${PWA2000_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# executables
make_executable(testParticle         testParticle.cc         "${AMPLITUDE_LIBS}")
make_executable(testDecayTopology    testDecayTopology.cc    "${AMPLITUDE_LIBS}")
make_executable(testWaveDescription  testWaveDescription.cc  "${AMPLITUDE_LIBS}")
make_executable(testWaveSet          testWaveSet.cc          "${AMPLITUDE_LIBS}")
make_executable(testWaveSetGenerator testWaveSetGenerator.cc "${AMPLITUDE_LIBS}")
make_executable(testAmplitude        testAmplitude.cc        "${AMPLITUDE_LIBS}")
make_executable(testAmplitudeTree    testAmplitudeTree.cc    "${AMPLITUDE_LIBS}")
make_executable(testIntegral         testIntegral.cc         "${AMPLITUDE_LIBS}")
