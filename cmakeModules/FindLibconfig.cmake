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
#//      cmake module for finding Libconfig installation
#//      Libconfig installation location is defined by environment variable $LIBCONFIG
#//
#//      following variables are defined:
#//      Libconfig_VERSION          - Libconfig version
#//      Libconfig_MAJOR_VERSION    - Libconfig major version
#//      Libconfig_MINOR_VERSION    - Libconfig minor version
#//      Libconfig_SUBMINOR_VERSION - Libconfig patch level
#//      Libconfig_DIR              - Libconfig installation directory
#//      Libconfig_INCLUDE_DIR      - Libconfig header directory
#//      Libconfig_LIBRARY_DIR      - Libconfig library directory
#//      Libconfig_LIBS             - Libconfig library files
#//
#//      Example usage:
#//          find_package(Libconfig 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(Libconfig_FOUND        TRUE)
set(Libconfig_ERROR_REASON "")

set(Libconfig_VERSION)
set(Libconfig_MAJOR_VERSION)
set(Libconfig_MINOR_VERSION)
set(Libconfig_SUBMINOR_VERSION)
set(Libconfig_DIR)
set(Libconfig_INCLUDE_DIR)
set(Libconfig_LIBRARY_DIR)
set(Libconfig_LIBS)


# try to get the environment variable pointing to the Libconfig installation
# directory
set(Libconfig_DIR $ENV{LIBCONFIG})


# find the library
set(_Libconfig_LIBRARY_NAMES "config++")
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_library(Libconfig_LIBS
		NAMES ${_Libconfig_LIBRARY_NAMES}
		PATHS ${Libconfig_DIR}/lib
		NO_DEFAULT_PATH
		)
else()
	# search system-wide
	find_library(Libconfig_LIBS
		NAMES ${_Libconfig_LIBRARY_NAMES}
		)
endif()
if(Libconfig_LIBS)
	get_filename_component(Libconfig_LIBRARY_DIR ${Libconfig_LIBS} DIRECTORY)
else()
	set(Libconfig_FOUND FALSE)
	if(Libconfig_DIR)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig library '${_Libconfig_LIBRARY_NAMES}' "
			"in directory '${Libconfig_DIR}/lib'.")
	else()
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig library '${_Libconfig_LIBRARY_NAMES}' "
			"in any standard library directory.")
	endif()
endif()
unset(_Libconfig_LIBRARY_NAMES)


# find the include directory
set(_Libconfig_HEADER_FILE_NAME "libconfig.h++")
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_path(Libconfig_INCLUDE_DIR
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		PATHS ${Libconfig_DIR}/include
		NO_DEFAULT_PATH
		)
else()
	# search system-wide
	find_path(Libconfig_INCLUDE_DIR
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		)
endif()
if(NOT Libconfig_INCLUDE_DIR)
	set(Libconfig_FOUND FALSE)
	if(Libconfig_DIR)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig include file '${_Libconfig_HEADER_FILE_NAME}' "
			"in directory '${Libconfig_DIR}/include'.")
	else()
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig include file '${_Libconfig_HEADER_FILE_NAME}' "
			"in any standard include directory.")
	endif()
endif()


# try to get the version from the header file
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_file(_Libconfig_HEADER_FILE
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		PATHS ${Libconfig_DIR}/include
		NO_DEFAULT_PATH
		)
else()
	# search system-wide
	find_file(_Libconfig_HEADER_FILE
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		)
endif()
if(_Libconfig_HEADER_FILE)

	# parse version string
	file(STRINGS ${_Libconfig_HEADER_FILE} _Libconfig_VERSIONS
		REGEX "LIBCONFIGXX_VER_(MAJOR|MINOR|REVISION)")
	list(LENGTH _Libconfig_VERSIONS _NMB_Libconfig_VERSIONS)
	if(NOT _NMB_Libconfig_VERSIONS EQUAL 3)
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Cannot determine Libconfig version from file '${_Libconfig_HEADER_FILE}'.")
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
	set(Libconfig_VERSION
		"${Libconfig_MAJOR_VERSION}.${Libconfig_MINOR_VERSION}.${Libconfig_SUBMINOR_VERSION}")
	unset(_Libconfig_VERSIONS)
	unset(_NMB_Libconfig_VERSIONS)
endif()
unset(_Libconfig_HEADER_FILE)
unset(_Libconfig_HEADER_FILE_NAME)


# check the version
if(NOT Libconfig_VERSION)
	set(Libconfig_FOUND FALSE)
	set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Could not extract version for Libconfig.")
else()
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


# report result
if(Libconfig_FOUND)
	message(STATUS "Found Libconfig version ${Libconfig_VERSION} in '${Libconfig_DIR}'.")
	message(STATUS "Using Libconfig include directory '${Libconfig_INCLUDE_DIR}'.")
	message(STATUS "Using Libconfig library '${Libconfig_LIBS}'.")
else()
	if(Libconfig_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested Libconfig installation:${Libconfig_ERROR_REASON}")
	else()
		if(NOT Libconfig_FIND_QUIETLY)
			message(STATUS "Libconfig version ${Libconfig_FIND_VERSION}+ was not found:${Libconfig_ERROR_REASON}")
		endif()
	endif()
endif()


# make variables changeable
mark_as_advanced(
	Libconfig_VERSION
	Libconfig_MAJOR_VERSION
	Libconfig_MINOR_VERSION
	Libconfig_SUBMINOR_VERSION
	Libconfig_DIR
	Libconfig_INCLUDE_DIR
	Libconfig_LIBRARY_DIR
	Libconfig_LIBS
	)
