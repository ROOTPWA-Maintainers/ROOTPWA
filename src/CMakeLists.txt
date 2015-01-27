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
#//      build file for RootPwa library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'src' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	SYSTEM ${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${CUDA_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	fitResult.cc
	complexMatrix.cc
	pwaLikelihood.cc
	pwaPlotter.cc
	massDepFit.cc
	massDepFitComponents.cc
	massDepFitFsmd.cc
	massDepFitLikeli.cc
	massDepFitModel.cc
	parameterSpace.cc
	)


# ROOT dictionary
set(ROOTPWA_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaDict.cc)
set(ROOTPWA_DICT_INCLUDE_DIRS
	${ROOT_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
root_generate_dictionary(
	"${ROOTPWA_DICTIONARY}"
	"${ROOTPWA_DICT_INCLUDE_DIRS}"
	"fitResult.h;pwaPlotter.h;complexMatrix.h"
	"rootPwaLinkDef.h"
	)
set(SOURCES ${SOURCES} ${ROOTPWA_DICTIONARY})


# library
set(THIS_LIB "RootPwa")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOT_LIBS}"
	"${CUDA_LIBS}"
	"${AMPLITUDE_LIBS}"
	)


# executables
make_executable(pwafit            pwafit.cc            ${THIS_LIB} "${UTILITIES_LIBS}")
make_executable(pwaplot           pwaplots.cc          ${THIS_LIB})

make_executable(pwaMassFit        pwaMassFit.cc        ${THIS_LIB} "${UTILITIES_LIBS}" "${Libconfig_LIBS}")

make_executable(stripIntegralsFromFitResult  stripIntegralsFromFitResult.cc ${ROOTPWA_LIBS} ${UTILITIES_LIBS})
make_executable(addIntegralsToFitResult      addIntegralsToFitResult.cc     ${ROOTPWA_LIBS} ${UTILITIES_LIBS})
make_executable(compareFitResultsIntegrals   compareFitResultsIntegrals.cc  ${ROOTPWA_LIBS} ${UTILITIES_LIBS})

# build tests
add_subdirectory(test)
