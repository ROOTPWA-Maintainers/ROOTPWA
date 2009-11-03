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

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pwa2000/libpp
                    ${CMAKE_SOURCE_DIR}/pwa2000/gamp)


SET(SOURCES integral.cc
            keyParse.cc
	    keyfile.cc
	    keyScan.cc
	    wave.cc
	    #keyScan.ll
	    #keyParse.yy
	    Tgamp.cc
)
	    
SET(MYFLAGS "-g -Wall")

SET(CMAKE_CXX_FLAGS ${MYFLAGS})
SET(CMAKE_CXX_LDFLAGS -g)

ADD_LIBRARY(pwa SHARED ${SOURCES})

#add libpp
SET(ROOTPWA_LIB "${LIBRARY_OUTPUT_PATH}")
SET(PWA2000_LIBS "${ROOTPWA_LIB}/libpp.so")


ADD_EXECUTABLE(gamp gamp.cc)
TARGET_LINK_LIBRARIES(gamp pwa)
TARGET_LINK_LIBRARIES(gamp ${PWA2000_LIBS})

ADD_EXECUTABLE(vamp vamp.cc)
TARGET_LINK_LIBRARIES(vamp pwa)
TARGET_LINK_LIBRARIES(vamp ${PWA2000_LIBS})


ADD_EXECUTABLE(int int.cc)
TARGET_LINK_LIBRARIES(int pwa)
TARGET_LINK_LIBRARIES(int ${PWA2000_LIBS})
