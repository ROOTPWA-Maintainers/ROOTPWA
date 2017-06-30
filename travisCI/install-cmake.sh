#!/bin/bash
set -exv

cd ${TRAVIS_BUILD_DIR}/deps/

if [ -d cmake-3.8.2-Linux-x86_64 ] ; then
	echo "CMake installation found in 'cmake-3.8.2-Linux-x86_64', using that."
else
	echo "No CMake installation found, installing a fresh one."

	wget https://cmake.org/files/v3.8/cmake-3.8.2-Linux-x86_64.tar.gz
	tar -xzf cmake-3.8.2-Linux-x86_64.tar.gz
	rm -rf cmake-3.8.2-Linux-x86_64.tar.gz
fi
ln -sfn cmake-3.8.2-Linux-x86_64 cmake
