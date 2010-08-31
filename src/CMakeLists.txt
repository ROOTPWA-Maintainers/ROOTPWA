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
#//      build file for RootPwa library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# set include directories
set(INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
  ${PWA2000_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
  TFitBin.cc
  TFitResult.cc
  fitResult.cc
  TCMatrix.cc
  TPWALikelihood.h
  TPWALikelihoodC.cc
  TPWAAmp.cc
  TPWAPhaseConstraint.cc
  TLogMultiGaus.cc
  TCovEllipse.cc
  utilities.cc
  pwaPlotter.cc
)

	    
# ROOT dictionary
set(ROOTPWA_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaDict.cc)
root_generate_dictionary(
  "${ROOTPWA_DICTIONARY}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "TFitBin.h;TFitResult.h;fitResult.h;pwaPlotter.h;TCMatrix.h;TCovEllipse.h"
  "rootPwaLinkDef.h"
)
set(SOURCES ${SOURCES} ${ROOTPWA_DICTIONARY})


# library
set(THIS_LIB "RootPwa")
make_shared_library("${THIS_LIB}" "${SOURCES}"
  "${ROOT_LIBS}"
  "${PWA2000_LIBS}"
)


# executables
make_executable(listwaves         listwaves.cc         ${THIS_LIB})
make_executable(pwafit            pwafit.cc            ${THIS_LIB})
make_executable(pwaConstraintFit  pwaConstraintFit.cc  ${THIS_LIB})
make_executable(selector          selector.cc          ${THIS_LIB})
make_executable(pwaplot           pwaplots.cc          ${THIS_LIB})

make_executable(testPWALikeli     testPWALikeli.cc     ${THIS_LIB})
make_executable(testPWAConstraint testPWAConstraint.cc ${THIS_LIB})
make_executable(testPWAAmp        testPWAAmp.cc        ${THIS_LIB})
make_executable(testLikelihood    testLikelihood.cc    ${THIS_LIB})


# test executables
add_test(testBuild        ${CMAKE_BUILD_TOOL})
add_test(testPWAAmp       ${EXECUTABLE_OUTPUT_PATH}/testPWAAmp)
add_test(testPWALikeli    ${CMAKE_CURRENT_SOURCE_DIR}/testPWALikeli.sh)
add_test(testPWAConstrain ${CMAKE_CURRENT_SOURCE_DIR}/testPWAConstraint.sh)
# add_test(testPWAFit       ${CMAKE_CURRENT_SOURCE_DIR}/testPWAFit.sh)
