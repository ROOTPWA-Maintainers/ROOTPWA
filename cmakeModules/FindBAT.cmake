#//-------------------------------------------------------------------------
#//
#// Description:
#//      cmake module for finding BAT installation
#//      BAT installation location is defined by environment variable
#//      $BATINSTALLDIR. If this variable does not exist, bat-config is
#//      searched in the path.
#//
#//      following variables are defined:
#//      BAT_VERSION      - BAT version
#//      BAT_ROOT_DIR     - BAT installation directory
#//      BAT_INCLUDE_DIR  - BAT header directory
#//      BAT_LIBRARIES    - BAT libraries
#//      BAT_CXX_FLAGS    - extra compiler flags when using BAT
#//      BAT_LINKER_FLAGS - extra linker flags required to link against BAT
#//
#//      Example usage:
#//          find_package(BAT 0.9.3 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(BAT_FOUND        TRUE)
set(BAT_ERROR_REASON "")

set(BAT_VERSION      NOTFOUND)
set(BAT_ROOT_DIR     NOTFOUND)
set(BAT_INCLUDE_DIR  NOTFOUND)
set(BAT_LIBRARIES    NOTFOUND)
set(BAT_CXX_FLAGS    "")
set(BAT_LINKER_FLAGS "")


# try to get the environment variable pointing to the BAT installation
# directory
set(BAT_ROOT_DIR $ENV{BATINSTALLDIR})
# if the environment variable BATINSTALLDIR is empty, look whether 'bat-config'
# is somewhere in the path
if(NOT BAT_ROOT_DIR)
	find_program(_BAT_CONFIG_EXECUTABLE bat-config)
	if(_BAT_CONFIG_EXECUTABLE)
		execute_process(COMMAND ${_BAT_CONFIG_EXECUTABLE} --prefix
			OUTPUT_VARIABLE BAT_ROOT_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	endif()
endif()
if(NOT BAT_ROOT_DIR)
	set(BAT_FOUND FALSE)
	set(BAT_ERROR_REASON "${BAT_ERROR_REASON} cannot find BAT directory. Either environment variable BAT_ROOT_DIR is not set correctly or bat-config is not in path.")
endif()


if(BAT_ROOT_DIR)

	# find the library
	set(_BAT_LIBRARY_DIR "${BAT_ROOT_DIR}/lib")
	if(NOT EXISTS "${_BAT_LIBRARY_DIR}")
		set(BAT_FOUND FALSE)
		set(BAT_ERROR_REASON "${BAT_ERROR_REASON} Directory '${_BAT_LIBRARY_DIR}' does not exist.")
	else()
		set(_BAT_LIBRARY_NAME "BAT")
		find_library(BAT_LIBRARIES
			NAMES ${_BAT_LIBRARY_NAME}
			PATHS "${_BAT_LIBRARY_DIR}"
			NO_DEFAULT_PATH)
		if(NOT BAT_LIBRARIES)
			set(BAT_FOUND FALSE)
			set(BAT_ERROR_REASON "${BAT_ERROR_REASON} Cannot find BAT library '${_BAT_LIBRARY_NAME}' in '${_BAT_LIBRARY_DIR}'.")
		endif()
		unset(_BAT_LIBRARY_NAME)
	endif()
	unset(_BAT_LIBRARY_DIR)

	# find the include directory
	set(BAT_INCLUDE_DIR "${BAT_ROOT_DIR}/include")
	if(NOT EXISTS "${BAT_INCLUDE_DIR}")
		set(BAT_FOUND FALSE)
		set(BAT_ERROR_REASON "${BAT_ERROR_REASON} Directory '${BAT_INCLUDE_DIR}' does not exist.")
	endif()

	# find bat-config executable
	# if it was already found above, this will not change anything, if
	# it was not found only look in BAT_ROOT_DIR
	find_program(_BAT_CONFIG_EXECUTABLE bat-config
		PATHS "${BAT_ROOT_DIR}/bin"
		NO_DEFAULT_PATH)

endif()


if(_BAT_CONFIG_EXECUTABLE)
	# prefer to extract version and OpenMP info from 'bat-config'

	# get version
	execute_process(COMMAND ${_BAT_CONFIG_EXECUTABLE} --version
		OUTPUT_VARIABLE BAT_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	# check whether BAT needs OpenMP
	execute_process(COMMAND ${_BAT_CONFIG_EXECUTABLE} --cflags
		OUTPUT_VARIABLE _BAT_CFLAGS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	string(FIND ${_BAT_CFLAGS} "-fopenmp" _BAT_OPENMP)
	if(NOT _BAT_OPENMP EQUAL -1)
		set(BAT_CXX_FLAGS "${BAT_CXX_FLAGS} -fopenmp")
	endif()
	unset(_BAT_CFLAGS)
	execute_process(COMMAND ${_BAT_CONFIG_EXECUTABLE} --libs
		OUTPUT_VARIABLE _BAT_LDFLAGS
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	string(FIND ${_BAT_LDFLAGS} "-fopenmp" _BAT_OPENMP)
	if(NOT _BAT_OPENMP EQUAL -1)
		set(BAT_LINKER_FLAGS "${BAT_LINKER_FLAGS} -fopenmp")
	endif()
	unset(_BAT_LDFLAGS)
	unset(_BAT_OPENMP)

else()
# extract version and OpenMP info from config.h

	set(_BAT_CONFIG_HEADER_FILE_NAME "${BAT_ROOT_DIR}/config.h")
	if(NOT EXISTS "${_BAT_CONFIG_HEADER_FILE_NAME}")
		set(BAT_FOUND FALSE)
		set(BAT_ERROR_REASON "${BAT_ERROR_REASON} BAT configuration '${_BAT_CONFIG_HEADER_FILE_NAME}' does not exist.")
	else()
		parse_version_from_configh_file("${_BAT_CONFIG_HEADER_FILE_NAME}" BAT_VERSION)

		# check whether BAT needs OpenMP
		file(STRINGS ${_BAT_CONFIG_HEADER_FILE_NAME} _BAT_OPENMP
			REGEX " THREAD_PARALLELIZATION ")
		string(REGEX REPLACE
			"#define THREAD_PARALLELIZATION ([01])"
			"\\1" _BAT_OPENMP_FLAG "${_BAT_OPENMP}")
		if(_BAT_OPENMP_FLAG STREQUAL _BAT_OPENMP)
			set(BAT_FOUND FALSE)
			set(BAT_ERROR_REASON "${BAT_ERROR_REASON} Could not extract whether BAT was build with OpenMP support from string '${_BAT_OPENMP}'.")
		endif()
		if(_BAT_OPENMP_FLAG)
			set(BAT_CXX_FLAGS "${BAT_CXX_FLAGS} -fopenmp")
			set(BAT_LINKER_FLAGS "${BAT_LINKER_FLAGS} -fopenmp")
		endif()
		unset(_BAT_OPENMP)
		unset(_BAT_OPENMP_FLAG)
	endif()
	unset(_BAT_CONFIG_HEADER_FILE_NAME)

endif()

# remove leading and trailing whitespaces
string(STRIP "${BAT_CXX_FLAGS}" BAT_CXX_FLAGS)
string(STRIP "${BAT_LINKER_FLAGS}" BAT_LINKER_FLAGS)


if(BAT_ERROR_REASON AND NOT BAT_FIND_QUIETLY)
	message(STATUS "Problems while finding the requested BAT installation:${BAT_ERROR_REASON}")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BAT
	FOUND_VAR BAT_FOUND
	REQUIRED_VARS BAT_ROOT_DIR BAT_VERSION BAT_INCLUDE_DIR BAT_LIBRARIES
	VERSION_VAR BAT_VERSION)
# additional reporting
if(BAT_FOUND AND NOT BAT_FIND_QUIETLY)
	message(STATUS "Using BAT include directory '${BAT_INCLUDE_DIR}'.")
	message(STATUS "Using BAT libraries '${BAT_LIBRARIES}'.")
	message(STATUS "Using extra CXX_FLAGS for BAT '${BAT_CXX_FLAGS}'.")
	message(STATUS "Using extra LINKER_FLAGS for BAT '${BAT_LINKER_FLAGS}'.")
endif()


# hide variables from normal GUI
mark_as_advanced(
	BAT_VERSION
	BAT_ROOT_DIR
	BAT_INCLUDE_DIR
	BAT_LIBRARIES
	BAT_CXX_FLAGS
	BAT_LINKER_FLAGS
	)


if(NOT BAT_FOUND)
	unset(BAT_VERSION)
	unset(BAT_ROOT_DIR)
	unset(BAT_INCLUDE_DIR)
	unset(BAT_LIBRARIES)
	unset(BAT_CXX_FLAGS)
	unset(BAT_LINKER_FLAGS)
endif()
