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
#//      switches to Clang/LLVM toolchain
#//
#//      usage:
#//        cmake -DCMAKE_TOOLCHAIN_FILE=../cmakeModules/SwitchToLlvmToolChain.cmake ..
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(LLVM_BIN_DIR       "/opt/llvm/build/bin"    )
set(CMAKE_C_COMPILER   "${LLVM_BIN_DIR}/clang"  )
set(CMAKE_CXX_COMPILER "${LLVM_BIN_DIR}/clang++")


# switch off automatic compiler identification by cmake
#include(CMakeForceCompiler)
# 'Generic' removes -rdynamic from linker, which llvm-ld does not support
#set(CMAKE_SYSTEM_NAME Generic)
# define compiler executables
#CMAKE_FORCE_C_COMPILER  ("${LLVM_BIN_DIR}/clang"   "Clang")
#CMAKE_FORCE_CXX_COMPILER("${LLVM_BIN_DIR}/clang++" "Clang")

# see http://mhungerford.blogspot.com/2010/10/cmake-and-clangllvm-fun.html
# switch off automatic compiler identification by cmake
#include(CMakeForceCompiler)
# 'Generic' removes -rdynamic from linker, which llvm-ld does not support
#set(CMAKE_SYSTEM_NAME Generic)
# define compiler executables
#CMAKE_FORCE_C_COMPILER  ("${LLVM_BIN_DIR}/clang"   "Clang")
#CMAKE_FORCE_CXX_COMPILER("${LLVM_BIN_DIR}/clang++" "Clang")

#define other executables
# set(CMAKE_AR      "${LLVM_BIN_DIR}/llvm-ar"     CACHE INTERNAL STRING)
# set(CMAKE_LINKER  "${LLVM_BIN_DIR}/llvm-ld"     CACHE INTERNAL STRING)
# set(CMAKE_NM      "${LLVM_BIN_DIR}/llvm-nm"     CACHE INTERNAL STRING)
# set(CMAKE_OBJDUMP "${LLVM_BIN_DIR}/llvm-ojdump" CACHE INTERNAL STRING)
# set(CMAKE_RANLIB  "${LLVM_BIN_DIR}/llvm-ranlib" CACHE INTERNAL STRING)
# 
# 
# set(CMAKE_C_LINK_EXECUTABLE "${LLVM_BIN_DIR}/llvm-ld <OBJECTS> -o  <TARGET> <CMAKE_C_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES>")
# set(CMAKE_CXX_LINK_EXECUTABLE "${LLVM_BIN_DIR}/llvm-ld <OBJECTS> -o  <TARGET> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES>")
 
# set(CMAKE_FIND_ROOT_PATH "${LLVM_BIN_DIR}")
# set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
# set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

#set(CMAKE_C_LINKER "${LLVM_BIN_DIR}/llvm-ld")
#set(CMAKE_CXX_LINKER "${LLVM_BIN_DIR}/llvm-ld")


# http://stackoverflow.com/questions/7031126/switching-between-gcc-and-clang-llvm-using-cmake

# CMake honors the environment variables CC and CXX upon detecting the C
# and C++ compiler to use:
# 
# $ export CC=/usr/bin/clang
# $ export CXX=/usr/bin/clang++
# $ cmake ..
# -- The C compiler identification is Clang
# -- The CXX compiler identification is Clang
# 
# The compiler specific flags can be overridden by putting them into a
# system wide CMake file an pointing the CMAKE_USER_MAKE_RULES_OVERRIDE
# variable to it. Create a file ~/ClangOverrides.txt with the following
# contents:
# 
# SET (CMAKE_C_FLAGS_INIT                "-Wall -std=c99")
# SET (CMAKE_C_FLAGS_DEBUG_INIT          "-g")
# SET (CMAKE_C_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
# SET (CMAKE_C_FLAGS_RELEASE_INIT        "-O4 -DNDEBUG")
# SET (CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
# 
# SET (CMAKE_CXX_FLAGS_INIT                "-Wall")
# SET (CMAKE_CXX_FLAGS_DEBUG_INIT          "-g")
# SET (CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
# SET (CMAKE_CXX_FLAGS_RELEASE_INIT        "-O4 -DNDEBUG")
# SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
# 
# The suffix _INIT will make CMake initialize the corresponding
# *_FLAGS_* variable with the given value. Then invoke cmake in the
# following way:
# 
# $ cmake -DCMAKE_USER_MAKE_RULES_OVERRIDE=~/ClangOverrides.txt ..
# 
# Finally to force the use of the LLVM binutils, set the internal
# variable _CMAKE_TOOLCHAIN_PREFIX. This variable is honored by the
# CMakeFindBinUtils module:
# 
# $ cmake -D_CMAKE_TOOLCHAIN_PREFIX=llvm- ..
# 
# Putting this all together you can write a shell wrapper which sets up
# the environment variables CC and CXX and then invokes cmake with the
# mentioned variable overrides.


#include(CMakePrintSystemInformation)
