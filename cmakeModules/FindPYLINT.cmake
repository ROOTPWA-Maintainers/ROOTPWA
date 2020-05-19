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
#//      PYLINT_FOUND      - indicates whether pylint was found
#//      PYLINT_VERSION    - version of pylint executable
#//      PYLINT_EXECUTABLE - path to pylint executable
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


set(PYLINT_FOUND TRUE)
set(PYLINT_ERROR_REASON "")

set(PYLINT_VERSION    NOTFOUND)
set(PYLINT_EXECUTABLE NOTFOUND)


# find pylint executable
find_program(PYLINT_EXECUTABLE pylint)
if(NOT PYTHON_EXECUTABLE)
	set(PYLINT_FOUND FALSE)
	set(PYLINT_ERROR_REASON "${PYLINT_ERROR_REASON} Did not find 'pylint' executable.")
endif()

# get version
if(PYLINT_EXECUTABLE)
	execute_process(COMMAND ${PYLINT_EXECUTABLE} --version
		OUTPUT_VARIABLE _PYLINT_VERSION_RAW
		ERROR_QUIET)
	parse_version_from_multline_string("${_PYLINT_VERSION_RAW}" "pylint[ \t]+" PYLINT_VERSION)
endif()


if(PYLINT_ERROR_REASON AND NOT PYLINT_FIND_QUIETLY)
	message(STATUS "Problems while finding the requested pylint installation:${PYLINT_ERROR_REASON}")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYLINT
	FOUND_VAR PYLINT_FOUND
	REQUIRED_VARS PYLINT_EXECUTABLE PYLINT_VERSION
	VERSION_VAR PYLINT_VERSION)
# additional reporting
if(PYLINT_FOUND AND NOT PYLINT_FIND_QUIETLY)
	message(STATUS "Using pylint executable at '${PYLINT_EXECUTABLE}'")
endif()


# hide variables from normal GUI
mark_as_advanced(
	PYLINT_VERSION
	PYLINT_EXECUTABLE
	)


if(NOT PYLINT_FOUND)
	unset(PYLINT_VERSION)
	unset(PYLINT_EXECUTABLE)
endif()
