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
	set(NEW_FILE_LIST)
	foreach(_OLD_FILE ${IN_FILE_LIST})
		string(REGEX REPLACE "^(.*)${OLD_EXT}$" "\\1${NEW_EXT}" _NEW_FILE ${_OLD_FILE})
		set(NEW_FILE_LIST ${NEW_FILE_LIST} ${_NEW_FILE})
	endforeach()
	unset(_OLD_FILE)
	set(${OUT_FILE_LIST} ${NEW_FILE_LIST})
endfunction(switch_file_extension)


# adds standard shared library
# additional libraries that should be linked to can be given as optional arguments
function(make_shared_library LIB_NAME SOURCES)
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


macro(enforce_out_of_source_build)
	if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
		message(FATAL_ERROR "Building this project in the source directory is not allowed. "
			"Please remove CMakeCache.txt, create a build directory, and run cmake there, for example:
        rm CMakeCache.txt
        mkdir build && cd build
        cmake ..")
	endif()
endmacro(enforce_out_of_source_build)


# removes item from property list of enabled features and adds it to the list of disabled ones
function(disable_feature FEATURE_NAME)
	if(DEBUG_OUTPUT)
		message(STATUS "disable_feature was called with the following arguments:
        FEATURE_NAME = '${FEATURE_NAME}'")
	endif()
	set(${FEATURE_NAME}_FOUND FALSE)
	get_property(_FEATURES GLOBAL PROPERTY ENABLED_FEATURES)
	if(DEBUG_OUTPUT)
		message(STATUS "Removing feature '${FEATURE_NAME}' from ENABLED_FEATURES = '${_FEATURES}'")
	endif()
	list(LENGTH _FEATURES _FEATURES_SIZE)
	if(_FEATURES_SIZE GREATER 0)
		list(REMOVE_ITEM _FEATURES ${FEATURE_NAME})
		set_property(GLOBAL PROPERTY ENABLED_FEATURES ${_FEATURES})
		set_property(GLOBAL APPEND PROPERTY DISABLED_FEATURES ${FEATURE_NAME})
	endif()
	unset(_FEATURES)
	unset(_FEATURES_SIZE)
endfunction(disable_feature)
