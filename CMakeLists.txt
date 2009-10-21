
# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

project(ROOTPWA)
ENABLE_TESTING()

# Build types
set(CMAKE_BUILD_TYPE None)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Werror")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")

set(CMAKE_CXX_FLAGS_NONE "-O3 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Werror")

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
add_subdirectory (tools)
add_subdirectory (src)
add_subdirectory (generators)
 
 


