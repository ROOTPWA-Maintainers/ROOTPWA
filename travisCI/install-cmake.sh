#!/bin/bash
set -exv

echo ">>> Using CMake version ${CMAKE_VERSION}"

CMAKE_VERSION_SHORT=$(echo "${CMAKE_VERSION}" | cut -d "." -f 1-2 -)

cd ${TRAVIS_BUILD_DIR}/deps/

if [ -d cmake-${CMAKE_VERSION} ]
then
	echo "    Existing CMake installation found in 'cmake-${CMAKE_VERSION}', using that."
else
	echo "    No CMake installation found, installing a fresh one."

	# download binary tarball
	wget https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz -O cmake.tar.gz

	# extract tarball
	mkdir cmake-${CMAKE_VERSION}
	tar -xzf cmake.tar.gz -C cmake-${CMAKE_VERSION} --strip-components=1
	rm -rf cmake.tar.gz
fi

ln -sfn cmake-${CMAKE_VERSION} cmake
