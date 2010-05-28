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
#//      main build file for PWA2000 libraries
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(PWA2000_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/pwa2000/libpp
  ${CMAKE_SOURCE_DIR}/pwa2000/gamp
  PARENT_SCOPE
)
set(PWA2000_LIBS
  ${LIBRARY_OUTPUT_PATH}/libpp.so
  ${LIBRARY_OUTPUT_PATH}/libpwa.so
  PARENT_SCOPE
)
add_subdirectory(libpp)
add_subdirectory(gamp)
