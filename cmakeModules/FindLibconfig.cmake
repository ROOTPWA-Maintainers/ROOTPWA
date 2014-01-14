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
#//      cmake module for finding libconfig installation
#//      libconfig installation location is defined by environment variable $LIBCONFIG
#//	 
#//      following variables are defined:
#//      Libconfig_VERSION          - libconfig version
#//      Libconfig_MAJOR_VERSION    - libconfig major version
#//      Libconfig_MINOR_VERSION    - libconfig minor version
#//      Libconfig_SUBMINOR_VERSION - libconfig patch level
#//      Libconfig_DIR              - libconfig installation directory
#//      Libconfig_INCLUDE_DIR      - libconfig header directory
#//      Libconfig_LIBRARY_DIR      - libconfig library directory
#//      Libconfig_LIBS             - libconfig library files
#//	 
#//      Example usage:
#//          find_package(Libconfig 1.4 REQUIRED)
#//	 
#//	 
#//-------------------------------------------------------------------------


set(Libconfig_FOUND        FALSE)
set(Libconfig_ERROR_REASON "")
set(Libconfig_DEFINITIONS  "")
set(Libconfig_LIBS)


set(Libconfig_DIR $ENV{LIBCONFIG})
if(NOT Libconfig_DIR)
	set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Environment variable LIBCONFIG='${Libconfig_DIR}' is not set correctly.")
else()

	set(Libconfig_FOUND TRUE)

	set(Libconfig_INCLUDE_DIR "${Libconfig_DIR}/include")
	if(NOT EXISTS "${Libconfig_INCLUDE_DIR}")
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Directory '${Libconfig_INCLUDE_DIR}' does not exist.")
	endif()

	set(Libconfig_LIBRARY_DIR "${Libconfig_DIR}/lib")
	if(NOT EXISTS "${Libconfig_LIBRARY_DIR}")
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Directory '${Libconfig_LIBRARY_DIR}' does not exist.")
	endif()
	
	set(_Libconfig_LIB_NAMES "config++")
	find_library(Libconfig_LIBS
		NAMES ${_Libconfig_LIB_NAMES}
		PATHS ${Libconfig_LIBRARY_DIR}
		NO_DEFAULT_PATH)
	if(NOT Libconfig_LIBS)
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Cannot find libconfig library '${_Libconfig_LIB_NAMES}' in '${Libconfig_LIBRARY_DIR}'.")
	endif()
	unset(_Libconfig_LIB_NAMES)

	set(_Libconfig_HEADER_FILE_NAME "libconfig.h++")
	find_file(_Libconfig_HEADER_FILE
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		PATHS ${Libconfig_INCLUDE_DIR}
		NO_DEFAULT_PATH)
	if(NOT _Libconfig_HEADER_FILE)
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Cannot find libconfig header file '${_Libconfig_HEADER_FILE_NAME}' in '${Libconfig_INCLUDE_DIR}'.")
	else()
		# parse version string
		file(STRINGS ${_Libconfig_HEADER_FILE} _Libconfig_VERSIONS
			REGEX "LIBCONFIGXX_VER_(MAJOR|MINOR|REVISION)")
		list(LENGTH _Libconfig_VERSIONS _NMB_Libconfig_VERSIONS)
		if(NOT _NMB_Libconfig_VERSIONS EQUAL 3)
			set(Libconfig_FOUND FALSE)
			set(Libconfig_ERROR_REASON "Cannot determine libconfig version from file '${_Libconfig_HEADER_FILE}'.")
		else()
			string(REGEX REPLACE
				"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_MAJOR[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
				"\\1"	Libconfig_MAJOR_VERSION "${_Libconfig_VERSIONS}")
			string(REGEX REPLACE
				"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_MINOR[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
				"\\1" Libconfig_MINOR_VERSION "${_Libconfig_VERSIONS}")
			string(REGEX REPLACE
				"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_REVISION[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
				"\\1"	Libconfig_SUBMINOR_VERSION "${_Libconfig_VERSIONS}")
		endif()
		unset(_Libconfig_VERSIONS)
		set(Libconfig_VERSION
			"${Libconfig_MAJOR_VERSION}.${Libconfig_MINOR_VERSION}.${Libconfig_SUBMINOR_VERSION}")
	endif()
	unset(_Libconfig_HEADER_FILE_NAME)
	unset(_Libconfig_HEADER_FILE)

	# compare version
	if(Libconfig_FIND_VERSION_EXACT)
		if(NOT Libconfig_VERSION VERSION_EQUAL Libconfig_FIND_VERSION)
			set(Libconfig_FOUND FALSE)
			set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Libconfig version ${Libconfig_VERSION} does not match requested version ${Libconfig_FIND_VERSION}.")
		endif()
	else()
		if(Libconfig_VERSION VERSION_LESS Libconfig_FIND_VERSION)
			set(Libconfig_FOUND FALSE)
			set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Libconfig version ${Libconfig_VERSION} is lower than requested version ${Libconfig_FIND_VERSION}.")
		endif()
	endif()

endif()


# make variables changeable
mark_as_advanced(
	Libconfig_INCLUDE_DIR
	Libconfig_LIBRARY_DIR
	Libconfig_LIBS
	Libconfig_DEFINITIONS
	)


# report result
if(Libconfig_FOUND)
	message(STATUS "Found libconfig version ${Libconfig_VERSION} in '${Libconfig_DIR}'.")
	message(STATUS "Using libconfig include directory '${Libconfig_INCLUDE_DIR}'.")
	message(STATUS "Using libconfig library '${Libconfig_LIBS}'.")
else()
	if(Libconfig_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested libconfig installation:${Libconfig_ERROR_REASON}")
	else()
		if(NOT Libconfig_FIND_QUIETLY)
			message(STATUS "libconfig version ${Libconfig_FIND_VERSION}+ was not found:${Libconfig_ERROR_REASON}")
		endif()
	endif()
endif()
