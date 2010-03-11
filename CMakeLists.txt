
# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

PROJECT(rootpwa)
ENABLE_TESTING()

# Build types
SET(CMAKE_BUILD_TYPE None)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Werror")
SET(CMAKE_CXX_LDFLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_NONE "-O3 -Wall -Werror")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Werror")

SET(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# Get ROOT Variables
SET(ROOT_CONFIG_SEARCHPATH $ENV{ROOTSYS}/bin)

FIND_PROGRAM(ROOT_CONFIG_EXECUTABLE NAMES root-config PATHS
 	   ${ROOT_CONFIG_SEARCHPATH}
 	   NO_DEFAULT_PATH)
	    	   
	IF(${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
 	  MESSAGE(FATAL_ERROR "ROOT not installed in the searchpath and ROOTSYS is not set. Please set ROOTSYS or add the path to your ROOT installation in the Macro FindROOT.cmake in th subdirectory cmake/modules.")
 	ELSE (${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND")
 	  STRING(REGEX REPLACE "(^.*)/bin/root-config" "\\1" test ${ROOT_CONFIG_EXECUTABLE})
 	  SET( ENV{ROOTSYS} ${test})
 	ENDIF (${ROOT_CONFIG_EXECUTABLE} MATCHES "ROOT_CONFIG_EXECUTABLE-NOTFOUND") 
 	
EXEC_PROGRAM( ${ROOT_CONFIG_EXECUTABLE}
 	        ARGS "--incdir"
 	        OUTPUT_VARIABLE root_headers )
 	      SET(ROOT_INCLUDE_DIR ${root_headers})
 	      # CACHE INTERNAL "")
 		
EXEC_PROGRAM( ${ROOT_CONFIG_EXECUTABLE}
 	        ARGS "--libdir"
 	        OUTPUT_VARIABLE root_libs )
 	      SET(ROOT_LIB ${root_libs})
 	      # CACHE INTERNAL "")



# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
add_subdirectory(pwa2000)
add_subdirectory(tools)
add_subdirectory(src)
add_subdirectory(generators)
add_subdirectory(keygen)
add_subdirectory(amplitude)
#add_subdirectory(gui)
 
