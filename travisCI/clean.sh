#!/bin/bash
set -exv

echo ">>> Cleaning ${TRAVIS_BUILD_DIR}/deps/"

cd "${TRAVIS_BUILD_DIR}"/deps/
ls -la .

# remove remainders of CMake installation (downloaded file)
rm -rf cmake.tar.gz

# remove remainders of ROOT installation (downloaded file)
rm -rf root.tar.gz
