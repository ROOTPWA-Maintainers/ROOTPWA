#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2016
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
#//
#// Description:
#//      cmake module for finding pylint executable
#//
#//      following variables are defined:
#//      PYLINT_FOUND              - pylint found
#//      PYLINT_EXECUTABLE         - path to pylint executable
#//      PYLINT_VERSION            - version of pylint executable
#//
#//      Example usage:
#//          find_package(PYLINT 1.5.5 Optional)
#//
#//
#// Author List:
#//      Sebastian Uhl        TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


unset(PYLINT_FOUND)
unset(PYLINT_EXECUTABLE)
unset(PYLINT_VERSION)

# find pylint executable
find_program(PYLINT_EXECUTABLE pylint)

# if pylint was found extract its version
if(PYLINT_EXECUTABLE)
	execute_process(COMMAND ${PYLINT_EXECUTABLE} --version
	                OUTPUT_VARIABLE _PYLINT_VERSION_RAW
	                ERROR_QUIET
	               )
	string(REGEX REPLACE "pylint ([0123456789\\.]+).*" "\\1" PYLINT_VERSION "${_PYLINT_VERSION_RAW}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYLINT
                                  REQUIRED_VARS PYLINT_EXECUTABLE PYLINT_VERSION
                                  VERSION_VAR PYLINT_VERSION
                                 )
