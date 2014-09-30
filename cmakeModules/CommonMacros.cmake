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
#//      collection of cmake macros
#//
#//
#// Author List:
#//      Boris Grube    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# takes list of file names and returns file name list with new extension
# example:
#   switch_file_extension("${CC_LIST}" ".cc" ".h" H_LIST)
function(switch_file_extension IN_FILE_LIST OLD_EXT NEW_EXT OUT_FILE_LIST)
	if(DEBUG_OUTPUT)
		message(STATUS "switch_file_extension was called with the following arguments:
        IN_FILE_LIST  = '${IN_FILE_LIST}'
        OLD_EXT       = '${OLD_EXT}'
        NEW_EXT       = '${NEW_EXT}'
        OUT_FILE_LIST = '${OUT_FILE_LIST}'")
	endif()
	set(_NEW_FILE_LIST)
	foreach(_OLD_FILE ${IN_FILE_LIST})
		string(REGEX REPLACE "^(.*)${OLD_EXT}$" "\\1${NEW_EXT}" _NEW_FILE ${_OLD_FILE})
		set(_NEW_FILE_LIST ${_NEW_FILE_LIST} ${_NEW_FILE})
	endforeach()
	set(${OUT_FILE_LIST} "${_NEW_FILE_LIST}"  PARENT_SCOPE)
	unset(_OLD_FILE)
	unset(_NEW_FILE_LIST)
endfunction(switch_file_extension)


# adds standard shared library
# additional libraries that should be linked to can be given as optional arguments
function(make_shared_library LIB_NAME SOURCES)
	message(STATUS ">>> setting up shared library '${LIB_NAME}'")
	if(DEBUG_OUTPUT)
		message(STATUS "make_shared_library was called with the following arguments:
        LIB_NAME = '${LIB_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
	endif()
	add_library(${LIB_NAME} SHARED ${SOURCES})
	# proccess link libraries in additional arguments
	foreach(_LIB ${ARGN})
		target_link_libraries(${LIB_NAME} ${_LIB})
	endforeach()
	unset(_LIB)
endfunction(make_shared_library)


# adds standard executable
# additional libraries that should be linked to can be given as optional arguments
function(make_executable EXE_NAME SOURCES)
	message(STATUS ">>> setting up executable '${EXE_NAME}'")
	if(DEBUG_OUTPUT)
		message(STATUS "make_executable was called with the following arguments:
        EXE_NAME = '${EXE_NAME}'
        SOURCES  = '${SOURCES}'
        ARGN     = '${ARGN}'")
	endif()
	add_executable(${EXE_NAME} ${SOURCES})
	# proccess link libraries in additional arguments
	foreach(_LIB ${ARGN})
		target_link_libraries(${EXE_NAME} ${_LIB})
	endforeach()
	unset(_LIB)
endfunction(make_executable)


# protects against building project in source directory
macro(enforce_out_of_source_build)
	if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
		message(FATAL_ERROR "Building this project in the source directory is not allowed. "
			"Please remove CMakeCache.txt, create a build directory, and run cmake there, for example:
        rm CMakeCache.txt
        mkdir build && cd build
        cmake ..")
	endif()
endmacro(enforce_out_of_source_build)


# little helper macro that is used to print log message when entering subdirectories
macro(message_setup_this_dir)
	message(STATUS "")
	get_filename_component(_THIS_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
	message(STATUS ">>> Setting up '${_THIS_DIR}' directory.")
	unset(_THIS_DIR)
endmacro(message_setup_this_dir)


# little helper function that converts source directory name into shared-library name
function(lib_name_from_dir_name SRC_DIR LIB_NAME)
	get_filename_component(_SRC_DIR_NAME ${SRC_DIR} NAME)
	# get fist letter and make it upper case
	string(SUBSTRING "${_SRC_DIR_NAME}" 0 1 _FIST_LETTER)
	string(TOUPPER "${_FIST_LETTER}" _FIST_LETTER)
	# get name w/o first letter
	string(LENGTH ${_SRC_DIR_NAME} _SRC_DIR_NAME_LENGTH)
	math(EXPR _SRC_DIR_NAME_LENGTH ${_SRC_DIR_NAME_LENGTH}-1)
	string(SUBSTRING "${_SRC_DIR_NAME}" 1 ${_SRC_DIR_NAME_LENGTH} _REST)
	# construct library name
	set(${LIB_NAME} "RootPwa${_FIST_LETTER}${_REST}" PARENT_SCOPE)
	unset(_SRC_DIR_NAME)
	unset(_FIST_LETTER)
	unset(_REST)
endfunction(lib_name_from_dir_name)


# little helper function that returns library name for current source directory
function(lib_name_for_current_dir LIB_NAME)
	get_filename_component(_THIS_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
	lib_name_from_dir_name(${_THIS_DIR} _LIB_NAME)
	set(${LIB_NAME} "${_LIB_NAME}" PARENT_SCOPE)
	unset(_THIS_DIR)
	unset(_LIB_NAME)
endfunction(lib_name_for_current_dir)


# little helper function that converts source directory name into an upper-case variable-name component
function(var_name_from_dir_name SRC_DIR VAR_NAME)
	get_filename_component(_SRC_DIR_NAME ${SRC_DIR} NAME)
	string(TOUPPER "${_SRC_DIR_NAME}" _SRC_DIR_NAME)
	set(${VAR_NAME} "RPWA_${_SRC_DIR_NAME}" PARENT_SCOPE)
	unset(_SRC_DIR_NAME)
endfunction(var_name_from_dir_name)
