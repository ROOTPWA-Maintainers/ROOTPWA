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
#//      build file for amplitude library 
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# define project
PROJECT(rootpwa)
ENABLE_TESTING()


ADD_DEFINITIONS(-Wno-deprecated)  # suppress deprecated header warnings from boost libraries when using gcc 4.4


SET(INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/pwa2000/libpp
  ${CMAKE_SOURCE_DIR}/pwa2000/gamp
  ${CMAKE_SOURCE_DIR}/src
)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES})


# list of source files that are compiled into library
SET(SOURCES
  particleProperties.cc
  particleDataTable.cc
  particle.cc
  angMomCoupl.cc
  interactionVertex.cc
  interactionVertex2.cc
  fsVertex.cc
  diffractiveDissVertex.cc
  diffractiveDissVertex2.cc
  isobarDecayVertex.cc
  isobarDecayVertex2.cc
  decayTopology.cc
  decayTopology2.cc
  isobarDecayTopology.cc
  isobarHelicityAmplitude.cc
)

# define ROOT libraries
SET(ROOT_BIN $ENV{ROOTSYS}/bin)
SET(ROOT_LIBS
  ${ROOT_LIB}/libCore.so 
  ${ROOT_LIB}/libCint.so
  ${ROOT_LIB}/libRIO.so		
  ${ROOT_LIB}/libNet.so
  ${ROOT_LIB}/libHist.so
  ${ROOT_LIB}/libGraf.so
  ${ROOT_LIB}/libGraf3d.so 
  ${ROOT_LIB}/libGpad.so
  ${ROOT_LIB}/libTree.so
  ${ROOT_LIB}/libThread.so
  ${ROOT_LIB}/libRint.so
  ${ROOT_LIB}/libPostscript.so
  ${ROOT_LIB}/libMatrix.so 
  ${ROOT_LIB}/libPhysics.so
  ${ROOT_LIB}/libMathCore.so
  ${ROOT_LIB}/libMinuit2.so
  libdl.so
  libm.so
)

SET(PWA2000_LIBS
  ${LIBRARY_OUTPUT_PATH}/libpp.so
  ${LIBRARY_OUTPUT_PATH}/libpwa.so
)


#-------------------------------------------------
#------- The library -----------------------------
#-------------------------------------------------

SET(THIS_LIB RootPwaAmplitude)
ADD_LIBRARY(${THIS_LIB} SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(${THIS_LIB} ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(${THIS_LIB} ${BOOST_LIBS})
TARGET_LINK_LIBRARIES(${THIS_LIB} ${PWA2000_LIBS})


#-------------------------------------------------
#------- Main executables ------------------------
#-------------------------------------------------

SET(THIS_EXE testParticle)
ADD_EXECUTABLE(${THIS_EXE} ${THIS_EXE}.cc)
TARGET_LINK_LIBRARIES(${THIS_EXE} ${THIS_LIB})

SET(THIS_EXE testDecayTopology)
ADD_EXECUTABLE(${THIS_EXE} ${THIS_EXE}.cc)
TARGET_LINK_LIBRARIES(${THIS_EXE} ${THIS_LIB})

SET(THIS_EXE testDecayCombinatorics)
ADD_EXECUTABLE(${THIS_EXE} ${THIS_EXE}.cc)
TARGET_LINK_LIBRARIES(${THIS_EXE} ${THIS_LIB})

SET(THIS_EXE testAmplitude)
ADD_EXECUTABLE(${THIS_EXE} ${THIS_EXE}.cc)
TARGET_LINK_LIBRARIES(${THIS_EXE} ${THIS_LIB})


#-------------------------------------------------
#------- Install Commands ------------------------
#-------------------------------------------------

INSTALL(FILES librootPwaAmplitude.so DESTINATION lib)
