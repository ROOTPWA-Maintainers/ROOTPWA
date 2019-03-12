#!/bin/bash
set -exv

echo ">>> Building ROOTPWA"

cd "${TRAVIS_BUILD_DIR}"/build

# configure the build process
cmake "${TRAVIS_BUILD_DIR}"

# run 'make'
make --jobs=${NMB_JOBS}

echo ">>> Running 'make test'"
# - if 'make' was successful, run 'make test'
# - if 'make test' was not successful, print the log and fail the build
make test || (cat Testing/Temporary/LastTest.log && false)

echo ">>> Running 'testMC.sh'"
# - if 'make test' was successful, run 'testMC.sh'
cd "${TRAVIS_BUILD_DIR}"
# set some environment variables
export ROOTPWA=$(pwd -P)
export LD_LIBRARY_PATH=${ROOTPWA}/build/lib:${LD_LIBRARY_PATH}
export PATH=${ROOTPWA}/build/bin:${PATH}
export PYTHONPATH=${ROOTPWA}/build/pyLib:${PYTHONPATH}
# run the test itself
cd "${TRAVIS_BUILD_DIR}"/test
./testMC.sh
