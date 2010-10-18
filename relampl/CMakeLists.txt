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
#// $Rev:: 331                         $: revision of last commit
#// $Author:: bgrube                   $: author of last commit
#// $Date:: 2010-07-27 18:51:06 +0200 #$: date of last commit
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


add_definitions(-Wno-deprecated)  # suppress deprecated header warnings from boost libraries when using gcc 4.4

#set compiler flags
#set(CMAKE_CXX_FLAGS "-Wall")

# set include directories
set(INCLUDE_DIRECTORIES
 ${CMAKE_CURRENT_SOURCE_DIR}
)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
TFracNum.cc
ClebschGordanBox.cc
TJwfTensor.cc
TSpinWaveFunction.cc
TLSAmpl.cc
TFhh.cc
)


# library
set(THIS_LIB "RelAmpl")
make_shared_library("${THIS_LIB}" "${SOURCES}"
)


# executables
make_executable(calcAmpl CalcAmpl.cc ${THIS_LIB})

