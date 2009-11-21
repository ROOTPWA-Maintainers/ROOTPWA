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
  $ENV{ROOTSYS}/include
  $ENV{LIBCONFIG}/include
  ${CMAKE_SOURCE_DIR}/pwa2000/libpp
  ${CMAKE_SOURCE_DIR}/pwa2000/gamp
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/generators
)

include_directories(${INCLUDE_DIRECTORIES})

SET(SOURCES  TProductionAmp.cc
             TBWProductionAmp.cc
             TPWWeight.cc
             TDiffractivePhaseSpace.cc)
   
#add ROOT lib
SET(ROOT_LIB $ENV{ROOTSYS}/lib)
SET(ROOT_BIN $ENV{ROOTSYS}/bin)
SET(ROOT_LIBS ${ROOT_LIB}/libCore.so
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
              ${ROOT_LIB}/libMinuit.so
              libdl.so
              libm.so)
#              ${ROOT_LIB}/libMinuit2.so
#add PWA lib
SET(PWA_LIB ${LIBRARY_OUTPUT_PATH})
SET(PWA_LIBS ${PWA_LIB}/libpp.so
             ${PWA_LIB}/libpwa.so)

#add rootpwalib
SET(ROOTPWA_LIB ${LIBRARY_OUTPUT_PATH})
SET(ROOTPWA_LIBS ${ROOTPWA_LIB}/librootpwa.so)

SET(LIBCONFIG_LIB $ENV{LIBCONFIG}/lib)
SET(LIBCONFIG ${LIBCONFIG_LIB}/libconfig++.so)


#-------------------------------------------------
#------ Root Dictionary --------------------------
#-------------------------------------------------

FIND_PROGRAM(ROOT_CINT_EXECUTABLE
  NAMES rootcint
  PATHS ${ROOT_BIN}
  NO_DEFAULT_PATH
 )

ADD_CUSTOM_COMMAND(OUTPUT rootpwagenDict.h ${CMAKE_CURRENT_BINARY_DIR}/rootpwagenDict.cc
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f rootpwagenDict.cc -c -DHAVE_CONFIG_H -I${CMAKE_SOURCE_DIR}/generators rootpwatoolsLinkDef.h)


set(GEN_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootpwagenDict.cc)

#-------------------------------------------------
#------- The library ---- ------------------------
#-------------------------------------------------

SET(SOURCES ${SOURCES} ${TOOLS_DICTIONARY})
ADD_LIBRARY(rootpwagen SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(rootpwagen ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(rootpwagen ${PWA_LIBS})
TARGET_LINK_LIBRARIES(rootpwagen ${ROOTPWA_LIBS})

#-------------------------------------------------
#------- Main executables ------------------------
#-------------------------------------------------

ADD_EXECUTABLE(genpw genpw.cc)
TARGET_LINK_LIBRARIES(genpw ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(genpw ${PWA_LIBS})
TARGET_LINK_LIBRARIES(genpw ${LIBCONFIG})
TARGET_LINK_LIBRARIES(genpw rootpwagen)
TARGET_LINK_LIBRARIES(genpw ${ROOTPWA_LIBS})

ADD_EXECUTABLE(evtweight evtweight.cc)
TARGET_LINK_LIBRARIES(evtweight ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(evtweight ${PWA_LIBS})
TARGET_LINK_LIBRARIES(evtweight ${ROOTPWA_LIBS})

ADD_EXECUTABLE(plotBWAmp plotBWAmp.cc)
TARGET_LINK_LIBRARIES(plotBWAmp ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(plotBWAmp rootpwagen)
#-------------------------------------------------
#------- Install options -------------------------
#-------------------------------------------------

#INSTALL(FILES DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
