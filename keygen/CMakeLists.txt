#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2009 Sebastian Neubert
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


PROJECT(rootpwa)


ENABLE_TESTING()


SET(INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/pwa2000/libpp
  ${CMAKE_SOURCE_DIR}/pwa2000/gamp
  ${CMAKE_SOURCE_DIR}/src
)

INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES})

# list of source files that are compiled into library
SET(SOURCES
  particleKey.cc
  waveKey.cc
  genKeyHelper.cc
)

# define ROOT libraries
#SET(ROOT_LIB $ENV{ROOTSYS}/lib)
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

# define PWA libraries
SET(PWA_LIB ${LIBRARY_OUTPUT_PATH})
SET(PWA_LIBS
  ${PWA_LIB}/libpp.so
  ${PWA_LIB}/libpwa.so
)

# add rootpwalib
SET(ROOTPWA_LIB ${LIBRARY_OUTPUT_PATH})
SET(ROOTPWA_LIBS ${ROOTPWA_LIB}/librootpwa.so)


#-------------------------------------------------
#------- The library -----------------------------
#-------------------------------------------------

SET(SOURCES
  ${SOURCES}
)

ADD_LIBRARY(rootpwaKey SHARED ${SOURCES})

TARGET_LINK_LIBRARIES(rootpwaKey ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(rootpwaKey ${PWA_LIBS})
TARGET_LINK_LIBRARIES(rootpwaKey ${ROOTPWA_LIBS})


#-------------------------------------------------
#------- Install Commands ------------------------
#-------------------------------------------------

INSTALL(FILES librootpwaKey.so DESTINATION lib)
