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
#//
#// Description:
#//      cmake module for finding NLopt installation
#//      NLopt installation location is defined by environment variable $NLOPT
#//
#//      following variables are defined:
#//      NLopt_DIR              - NLopt installation directory
#//      NLopt_INCLUDE_DIR      - NLopt header directory
#//      NLopt_LIBRARY_DIR      - NLopt library directory
#//      NLopt_LIBS             - NLopt library files
#//
#//      Example usage:
#//          find_package(NLopt 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(NLopt_FOUND        FALSE)
set(NLopt_ERROR_REASON "")
set(NLopt_DEFINITIONS  "")
set(NLopt_LIBS)
set(NLopt_DIR $ENV{NLOPT})


if(NOT NLopt_DIR)
	# search in system directories if environment variable NLOPT is not set

	set(NLopt_FOUND TRUE)

	set(_NLopt_LIB_NAMES "nlopt")
	find_library(NLopt_LIBS
		NAMES ${_NLopt_LIB_NAMES})
	if(NOT NLopt_LIBS)
		set(NLopt_FOUND FALSE)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Cannot find NLopt library '${_NLopt_LIB_NAMES}'.")
	else()
		get_filename_component(NLopt_DIR ${NLopt_LIBS} PATH)
	endif()
	unset(_NLopt_LIB_NAMES)

	set(_NLopt_HEADER_FILE_NAME "nlopt.hpp")
	find_file(_NLopt_HEADER_FILE
		NAMES ${_NLopt_HEADER_FILE_NAME})
	if(NOT _NLopt_HEADER_FILE)
		set(NLopt_FOUND FALSE)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Cannot find NLopt header file '${_NLopt_HEADER_FILE_NAME}'.")
	endif()
	unset(_NLopt_HEADER_FILE_NAME)
	unset(_NLopt_HEADER_FILE)

	if(NOT NLopt_FOUND)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} NLopt not found in system directories (and environment variable NLOPT is not set).")
	endif()

else()
	# search in directory defined by environment variable NLOPT

	set(NLopt_FOUND TRUE)

	set(NLopt_INCLUDE_DIR "${NLopt_DIR}/include")
	if(NOT EXISTS "${NLopt_INCLUDE_DIR}")
		# also check directory of CMake build introduced in 2.5.0
		set(NLopt_INCLUDE_DIR "${NLopt_DIR}/build/src/api")
		if(NOT EXISTS "${NLopt_INCLUDE_DIR}")
			set(NLopt_FOUND FALSE)
			set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Neither directory '${NLopt_DIR}/include' nor '${NLopt_DIR}/build/src/api' exists.")
		endif()
	endif()

	set(NLopt_LIBRARY_DIR "${NLopt_DIR}/lib")
	if(NOT EXISTS "${NLopt_LIBRARY_DIR}")
		# also check directory of CMake build introduced in 2.5.0
		set(NLopt_LIBRARY_DIR "${NLopt_DIR}/build")
		if(NOT EXISTS "${NLopt_LIBRARY_DIR}")
			set(NLopt_FOUND FALSE)
			set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Neither directory '${NLopt_DIR}/lib' nor '${NLopt_DIR}/build' exists.")
		endif()
	endif()

	set(_NLopt_LIB_NAMES "nlopt_cxx")
	find_library(NLopt_LIBS
		NAMES ${_NLopt_LIB_NAMES}
		PATHS ${NLopt_LIBRARY_DIR}
		NO_DEFAULT_PATH)
	if(NOT NLopt_LIBS)
		set(NLopt_FOUND FALSE)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Cannot find NLopt library '${_NLopt_LIB_NAMES}' in '${NLopt_LIBRARY_DIR}'.")
	endif()
	unset(_NLopt_LIB_NAMES)

	set(_NLopt_HEADER_FILE_NAME "nlopt.hpp")
	find_file(_NLopt_HEADER_FILE
		NAMES ${_NLopt_HEADER_FILE_NAME}
		PATHS ${NLopt_INCLUDE_DIR}
		NO_DEFAULT_PATH)
	if(NOT _NLopt_HEADER_FILE)
		set(NLopt_FOUND FALSE)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Cannot find NLopt header file '${_NLopt_HEADER_FILE_NAME}' in '${NLopt_INCLUDE_DIR}'.")
	endif()
	unset(_NLopt_HEADER_FILE_NAME)
	unset(_NLopt_HEADER_FILE)

endif()


# make variables changeable
mark_as_advanced(
	NLopt_INCLUDE_DIR
	NLopt_LIBRARY_DIR
	NLopt_LIBS
	NLopt_DEFINITIONS
	)


# report result
if(NLopt_FOUND)
	message(STATUS "Found NLopt in '${NLopt_DIR}'.")
	message(STATUS "Using NLopt include directory '${NLopt_INCLUDE_DIR}'.")
	message(STATUS "Using NLopt library '${NLopt_LIBS}'.")
else()
	if(NLopt_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested NLopt installation:${NLopt_ERROR_REASON}")
	else()
		if(NOT NLopt_FIND_QUIETLY)
			message(STATUS "NLopt was not found:${NLopt_ERROR_REASON}")
		endif()
	endif()
endif()
