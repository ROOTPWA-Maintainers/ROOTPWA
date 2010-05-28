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
#//      build file for RootPwaGen
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
  ${LIBCONFIG_INCLUDE_DIR}
  ${PWA2000_INCLUDE_DIR}
  ${ROOTPWA_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
  nBodyPhaseSpaceGen.cc
  TProductionAmp.cc
  TBWProductionAmp.cc
  TPWWeight.cc
  TDiffractivePhaseSpace.cc
)
   

# ROOT dictionary
set(GEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaGenDict.cc)
root_generate_dictionary(
  "${GEN_DICTIONARY}"
  "${CMAKE_CURRENT_SOURCE_DIR};${ROOTPWA_INCLUDE_DIR}"
  "nBodyPhaseSpaceGen.h"
  "rootPwaGenLinkDef.h"
)
set(SOURCES ${SOURCES} ${GEN_DICTIONARY})


# library
set(THIS_LIB "RootPwaGen")
make_shared_library("${THIS_LIB}" "${SOURCES}"
  "${ROOT_LIBS}"
  "${PWA2000_LIBS}"
  "${ROOTPWA_LIBS}"
)


# executables
make_executable(genpw     genpw.cc     ${THIS_LIB} "${LIBCONFIG_LIBS}")
make_executable(evtweight evtweight.cc ${THIS_LIB})
make_executable(plotBWAmp plotBWAmp.cc ${THIS_LIB})
