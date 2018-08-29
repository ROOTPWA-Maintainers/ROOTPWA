#!/bin/bash
set -exv

CMAKE_VERSION=3.12.1
CMAKE_VERSION_SHORT=$(echo "${CMAKE_VERSION}" | cut -d "." -f 1-2 -)

cd ${TRAVIS_BUILD_DIR}/deps/

if [ -d cmake-${CMAKE_VERSION}-Linux-x86_64 ] ; then
	echo "CMake installation found in 'cmake-${CMAKE_VERSION}-Linux-x86_64', using that."
else
	echo "No CMake installation found, installing a fresh one."

	wget https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
	tar -xzf cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
        rm -rf cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz
fi
ln -sfn cmake-${CMAKE_VERSION}-Linux-x86_64 cmake
