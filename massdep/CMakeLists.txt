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
#//      build file for RootPwaMassDep
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'massdep' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${Libconfig_INCLUDE_DIR}
	${ROOTPWATOOLS_INCLUDE_DIR}
	${GEN_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${MASSDEP_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	/home/sneubert/Projects/boost_1_45_0/boost/numeric/ublas/
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	decayChannels.cc
	#mcPhaseSpace.cc
	cauchyIntegral.cc
	dMatrixAmp.cc
	dMatrixSet.cc
	#dynMassDep.cc
	)


# ROOT dictionary
#set(GEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaMassDepDict.cc)
#root_generate_dictionary(
#  "${GEN_DICTIONARY}"
#  "${CMAKE_CURRENT_SOURCE_DIR};${ROOTPWA_INCLUDE_DIR}"
#  ""
#  "rootPwaMassDepLinkDef.h"
#)
#set(SOURCES ${SOURCES} ${GEN_DICTIONARY})


# library
set(THIS_LIB "RootPwaMassDep")
set(GSL_LIBS
	libgsl.so
	libgslcblas.so
	)
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${ROOTPWAGEN_LIBS}"
	"${GSL_LIBS}"
	)



# executables
#make_executable(testMassDep testMassDep.cc ${THIS_LIB} ${ROOT_LIBS} )
make_executable(testKMatrix testKMatrix.cc ${THIS_LIB} ${ROOT_LIBS} )
make_executable(testRescat  testRescat.cc  ${THIS_LIB} ${ROOT_LIBS} ${GSL_LIBS})
make_executable(testDMatrix testDMatrix.cc ${THIS_LIB} ${ROOT_LIBS} ${GSL_LIBS})
make_executable(mixingRes   mixingRes.cc   ${THIS_LIB} ${ROOT_LIBS} ${GSL_LIBS})
