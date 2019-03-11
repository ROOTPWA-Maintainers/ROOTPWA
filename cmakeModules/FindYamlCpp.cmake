#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2014
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
#//      cmake module for finding yaml-cpp installation
#//      yaml-cpp installation location is defined by environment variable
#//      $YAML_CPP
#//
#//      following variables are defined:
#//      YamlCpp_VERSION     - yaml-cpp version
#//      YamlCpp_DIR         - yaml-cpp installation directory
#//      YamlCpp_INCLUDE_DIR - yaml-cpp header directory
#//      YamlCpp_LIBRARY_DIR - yaml-cpp library directory
#//      YamlCpp_LIBS        - yaml-cpp library files
#//
#//      Example usage:
#//          find_package(YamlCpp 0.5 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(YamlCpp_FOUND        TRUE)
set(YamlCpp_ERROR_REASON "")

set(YamlCpp_VERSION     NOTFOUND)
set(YamlCpp_DIR         NOTFOUND)
set(YamlCpp_INCLUDE_DIR NOTFOUND)
Set(YamlCpp_LIBRARY_DIR NOTFOUND)
set(YamlCpp_LIBS        NOTFOUND)


# try to get the environment variable pointing to the yaml-cpp installation
# directory
set(YamlCpp_DIR $ENV{YAML_CPP})


# find the library
set(_YamlCpp_LIBRARY_NAME "yaml-cpp")
if(YamlCpp_DIR)
	# search only in YamlCpp_DIR
	find_library(YamlCpp_LIBS
		NAMES ${_YamlCpp_LIBRARY_NAME}
		PATHS "${YamlCpp_DIR}/lib"
		      "${YamlCpp_DIR}/build"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_library(YamlCpp_LIBS
		NAMES ${_YamlCpp_LIBRARY_NAME})
endif()
if(YamlCpp_LIBS)
	get_filename_component(YamlCpp_LIBRARY_DIR ${YamlCpp_LIBS} DIRECTORY)
else()
	set(YamlCpp_FOUND FALSE)
	if(YamlCpp_DIR)
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp library '${_YamlCpp_LIBRARY_NAME}' "
			"in directories '${YamlCpp_DIR}/lib' or '${YamlCpp_DIR}/build'.")
	else()
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp library '${_YamlCpp_LIBRARY_NAME}' "
			"in any standard library directory.")
	endif()
endif()
unset(_YamlCpp_LIBRARY_NAME)


# find the include directory
set(_YamlCpp_HEADER_FILE_NAME "yaml-cpp/yaml.h")
if(YamlCpp_DIR)
	# search only in YamlCpp_DIR
	find_path(YamlCpp_INCLUDE_DIR
		NAMES ${_YamlCpp_HEADER_FILE_NAME}
		PATHS "${YamlCpp_DIR}/include"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_path(YamlCpp_INCLUDE_DIR
		NAMES ${_YamlCpp_HEADER_FILE_NAME})
endif()
if(NOT YamlCpp_INCLUDE_DIR)
	set(YamlCpp_FOUND FALSE)
	if(YamlCpp_DIR)
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp include file '${_YamlCpp_HEADER_FILE_NAME}' "
			"in directory '${YamlCpp_DIR}/include'.")
	else()
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp include file '${_YamlCpp_HEADER_FILE_NAME}' "
			"in any standard include directory.")
	endif()
endif()
unset(_YamlCpp_HEADER_FILE_NAME)


# try to get the version from the pkgconfig file
set(_YamlCpp_PKG_CONFIG_FILE_NAME "yaml-cpp.pc")
if(YamlCpp_DIR)
	# search only in YamlCpp_DIR
	find_file(_YamlCpp_PC_FILE
		NAMES ${_YamlCpp_PKG_CONFIG_FILE_NAME}
		PATHS "${YamlCpp_DIR}/lib/pkgconfig"
		      "${YamlCpp_DIR}/build"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_file(_YamlCpp_PC_FILE
		NAMES ${_YamlCpp_PKG_CONFIG_FILE_NAME}
		PATHS "${YamlCpp_LIBRARY_DIR}/pkgconfig"
		      "${YamlCpp_LIBRARY_DIR}")
endif()
unset(_YamlCpp_PKG_CONFIG_FILE_NAME)
if(_YamlCpp_PC_FILE)
	parse_version_from_pkg_config_file("${_YamlCpp_PC_FILE}" YamlCpp_VERSION)
endif()
unset(_YamlCpp_PC_FILE)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YamlCpp
	FOUND_VAR YamlCpp_FOUND
	REQUIRED_VARS YamlCpp_VERSION YamlCpp_INCLUDE_DIR YamlCpp_LIBRARY_DIR YamlCpp_LIBS
	VERSION_VAR YamlCpp_VERSION)
# additional reporting
if(YamlCpp_FOUND)
	message(STATUS "Using yaml-cpp include directory '${YamlCpp_INCLUDE_DIR}'.")
	message(STATUS "Using yaml-cpp library '${YamlCpp_LIBS}'.")
else()
	message(STATUS "Unable to find requested YamlCpp installation:${YamlCpp_ERROR_REASON}")
endif()


# hide variables from normal GUI
mark_as_advanced(
	YamlCpp_VERSION
	YamlCpp_DIR
	YamlCpp_INCLUDE_DIR
	YamlCpp_LIBRARY_DIR
	YamlCpp_LIBS
	)


if(NOT YamlCpp_FOUND)
	unset(YamlCpp_VERSION)
	unset(YamlCpp_DIR)
	unset(YamlCpp_INCLUDE_DIR)
	unset(YamlCpp_LIBRARY_DIR)
	unset(YamlCpp_LIBS)
endif()
