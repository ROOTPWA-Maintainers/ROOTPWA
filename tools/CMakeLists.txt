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
#//      build file for RootPwaTools library
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
  TPDGEntry.cc
  TPDGDB.cc
  FSParticle.cc
  NParticleState.cc
  NParticleEvent.cc
)
   

# ROOT dictionary
set(TOOLS_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaToolsDict.cc)
root_generate_dictionary(
  "${TOOLS_DICTIONARY}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "TPDGEntry.h;TPDGDB.h"
  "rootPwaToolsLinkDef.h"
)
set(SOURCES ${SOURCES} ${TOOLS_DICTIONARY})


# library
set(THIS_LIB "RootPwaTools")
make_shared_library("${THIS_LIB}" "${SOURCES}"
  "${ROOT_LIBS}"
  "${PWA2000_LIBS}"
)


# executables
make_executable(addamp          addamp.cc          ${THIS_LIB})
make_executable(reflect         reflect.cc         ${THIS_LIB})
make_executable(plotGampMassDep plotGampMassDep.cc ${THIS_LIB})
make_executable(key2tex         key2tex.cc         ${THIS_LIB})
make_executable(addisokey       addisokey.cc       ${THIS_LIB})
make_executable(evt2tree        evt2tree.cc        ${THIS_LIB})
make_executable(binevt          binevt.cc          ${THIS_LIB})
make_executable(mutator         mutator.cc         ${THIS_LIB})
make_executable(getThresholds   getThresholds.cc   ${THIS_LIB})
