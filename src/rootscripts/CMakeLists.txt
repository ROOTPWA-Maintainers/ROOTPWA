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
#// $Rev:: 550                         $: revision of last commit
#// $Author:: sneubert                 $: author of last commit
#// $Date:: 2010-12-10 15:16:56 +0100 #$: date of last commit
#//
#// Description:
#//      build file for RootPwaScripts library
#//
#//
#// Author List:
#//      Prometeusz Jasinski  KPH-MZ         (original author)
#//      Sebastian Neubert    TUM            (reused from)
#//      
#//
#//
#//-------------------------------------------------------------------------

 # gSystem->Load("libRootPwa.so");
 # gSystem->AddIncludePath("-I$ROOTPWA/utilities");
 # gSystem->AddIncludePath("-I$ROOTPWA/src");

# set include directories
set(INCLUDE_DIRECTORIES
#  ${Boost_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
  ${UTILITIES_INCLUDE_DIR}
  ${ROOTPWA_INCLUDE_DIR}
# ${CUDA_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	loadFitResult.C
	convertTFitResultTree.C
	plotIntensity.C
	plotPhase.C
	plotCoherence.C
	plotAllIntensities.C
	plotSpinTotals.C
	# plot4.C
)

	    
# ROOT dictionary
#set(ROOTPWA_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootPwaDict.cc)
#root_generate_dictionary(
#  "${ROOTPWA_DICTIONARY}"
#  "${INCLUDE_DIRECTORIES}"
#  "TFitBin.h;TFitResult.h;fitResult.h;pwaPlotter.h;TCMatrix.h;TCovEllipse.h"
#  "rootPwaLinkDef.h"
#)
#set(SOURCES ${SOURCES} ${ROOTPWA_DICTIONARY})


# library
set(THIS_LIB "RootPwaScripts")
make_shared_library("${THIS_LIB}" "${SOURCES}"
  "${ROOT_LIBS}"
  "${ROOTPWA_LIBS}"
)
