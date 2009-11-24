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

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/pwa2000/libpp)

SET(SOURCES event.cc
            massDep.cc
	    particle.cc
	    pputil.cc
	    lorentz.cc
	    matrix.cc
	    particleData.cc
	    Vec.cc
)
	    
SET(MYFLAGS "-O2 -Wall -Werror")

SET(CMAKE_CXX_FLAGS ${MYFLAGS})
SET(CMAKE_LDCXX_FLAGS "-O2")

ADD_LIBRARY(pp SHARED ${SOURCES})

#INSTALL(FILES libpp.a DESTINATION lib)
#INSTALL(FILES event.h massDep.h particleData.h pputil.h lorentz.h matrix.h particle.h Vec.h DESTINATION include)