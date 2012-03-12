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
#// $Rev:: 266                         $: revision of last commit
#// $Author:: promme                   $: author of last commit
#// $Date:: 2010-06-17 10:35:01 +0200 #$: date of last commit
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


# this is just a quick-hack to get rid of the compiler warnings
# "cast to pointer from integer of different size [-Werror=int-to-pointer-cast]"
# that are triggered by casting IDs to pointers
# !!!this should be fixed in the code!!!
# unfortunately gcc 4.4.5 gives an error using -Wno-int-to-pointer-cast
if(CMAKE_COMPILER_IS_GNUCXX AND GCC_VERSION VERSION_GREATER "4.6.0")
	add_definitions(-Wno-int-to-pointer-cast)
endif()

# set include directories
set(INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
  ${Libconfig_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTPWA_DICTIONARY}
  ${SCRIPTS_INCLUDE_DIR}
#  ${PWA2000_INCLUDE_DIR}
  ${ROOTPWA_INCLUDE_DIR}
  ${UTILITIES_INCLUDE_DIR}
  ${GEN_INCLUDE_DIR}
)
include_directories(${INCLUDE_DIRECTORIES})

# source files that are compiled into library
set(SOURCES
  TrpwaMainFrame.cc
  TrpwaSessionManager.cc
  TrpwaWaveSelectFrame.cc
  TrpwaJobManager.cc
  TrpwaEventTreeHandler.cc
  TrpwaCommonTools.cc
  TrpwaPlotAmpsFrame.cc
  TrpwaFitOptionsFrame.cc
) 

# ROOT dictionary
set(ROOTPWAGUI_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaGuiDict.cc)
root_generate_dictionary(
  "${ROOTPWAGUI_DICTIONARY}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "TrpwaMainFrame.h;TrpwaWaveSelectFrame.h;TrpwaPlotAmpsFrame.h;TrpwaFitOptionsFrame.h"
  "rootPwaGuiLinkDef.h"
)
# message(" using rootpwa_dictionary "${ROOTPWA_DICTIONARY})
# setting roopwa dictionary does not work yet since I'm too dumb to use cmake
set(SOURCES ${SOURCES} ${ROOTPWAGUI_DICTIONARY} ${ROOTPWA_DICTIONARY})

# message(" my sources are "${SOURCES})  


# library
set(THIS_LIB "RootPwaGui")
make_shared_library("${THIS_LIB}" "${SOURCES}"
  "${ROOT_LIBS}"
   "${Libconfig_LIBS}"
   "${ROOTPWASCRIPTS_LIBS}"
#  "${PWA2000_LIBS}"
   "${ROOTPWA_LIBS}"
   "${ROOTPWAGEN_LIBS}"
)


# executables
make_executable(rootpwa_gui rootpwa_gui.cpp ${THIS_LIB})
make_executable(rootpwa_fitresult_inspector_gui rootpwa_fitresult_inspector_gui.cpp ${THIS_LIB})
