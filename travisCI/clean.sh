#!/bin/bash
set -exv

echo ">>> Cleaning ${TRAVIS_BUILD_DIR}/deps/"

cd ${TRAVIS_BUILD_DIR}/deps/

# remove remainders of CMake installation (downloaded file)
rm -rf cmake-3.8.2-Linux-x86_64.tar.gz
# use CMake installed by Travis, remove the version in out cache
rm -rf cmake
rm -rf cmake-3.8.2-Linux-x86_64

# remove remainders of ROOT installation (downloaded file)
rm -rf root.tar.gz
