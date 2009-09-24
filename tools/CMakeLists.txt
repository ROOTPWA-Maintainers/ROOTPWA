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
  $ENV{PWA2000}/include
  ${CMAKE_SOURCE_DIR}/tools
)

include_directories(${INCLUDE_DIRECTORIES})

SET(SOURCES TPDGEntry.cc
            TPDGDB.cc)
   
#SET(MYFLAGS "-g")

#SET(CMAKE_CXX_FLAGS ${MYFLAGS})
#SET(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} -g")

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
SET(PWA_LIB $ENV{PWA2000}/lib)
SET(PWA_LIBS ${PWA_LIB}/libpp.a
             ${PWA_LIB}/libpwa.a)

#add rootpwalib
SET(ROOTPWA_LIB $ENV{ROOTPWA}/lib)
SET(ROOTPWA_LIBS ${ROOTPWA_LIB}/librootpwa.so)


#-------------------------------------------------
#------ Root Dictionary --------------------------
#-------------------------------------------------

FIND_PROGRAM(ROOT_CINT_EXECUTABLE
  NAMES rootcint
  PATHS ${ROOT_BIN}
  NO_DEFAULT_PATH
 )

ADD_CUSTOM_COMMAND(OUTPUT rootpwatoolsDict.h ${CMAKE_CURRENT_BINARY_DIR}/rootpwatoolsDict.cc
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f rootpwatoolsDict.cc -c -DHAVE_CONFIG_H -I${CMAKE_SOURCE_DIR}/tools TPDGEntry.h TPDGDB.h rootpwatoolsLinkDef.h)


set(TOOLS_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/rootpwatoolsDict.cc)

#-------------------------------------------------
#------- The library ---- ------------------------
#-------------------------------------------------

SET(SOURCES ${SOURCES} ${TOOLS_DICTIONARY})
ADD_LIBRARY(rootpwatools SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(rootpwatools ${ROOT_LIBS})

#-------------------------------------------------
#------- Main executables ------------------------
#-------------------------------------------------

ADD_EXECUTABLE(addamp addamp.cc)

ADD_EXECUTABLE(reflect reflect.cc)
TARGET_LINK_LIBRARIES(reflect ${PWA_LIBS})

ADD_EXECUTABLE(key2tex key2tex.cc)
TARGET_LINK_LIBRARIES(key2tex ${ROOT_LIBS})

ADD_EXECUTABLE(addisokey addisokey.cc)
TARGET_LINK_LIBRARIES(addisokey ${ROOT_LIBS})


ADD_EXECUTABLE(mutator mutator.cc)
TARGET_LINK_LIBRARIES(mutator ${ROOT_LIBS})

ADD_EXECUTABLE(getThresholds getThresholds.cc)
TARGET_LINK_LIBRARIES(getThresholds ${ROOT_LIBS})
TARGET_LINK_LIBRARIES(getThresholds ${PWA_LIBS})


#-------------------------------------------------
#------- Install options -------------------------
#-------------------------------------------------

INSTALL(FILES addamp key2tex reflect mutator addisokey DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
