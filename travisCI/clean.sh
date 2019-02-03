#!/bin/bash
set -exv

echo ">>> Cleaning ${TRAVIS_BUILD_DIR}/deps/"

cd "${TRAVIS_BUILD_DIR}"/deps/
echo ">>> Before cleaning"
ls -la .

# remove everything that is not in the list of excluded files/directories
EXCLUDES=("bat"   "bat-${BAT_VERSION}"
          "boost" "boost-${BOOST_VERSION}"
          "cmake" "cmake-${CMAKE_VERSION}"
          "root"  "root-${ROOT_VERSION}"
          "yaml"  "yaml-${YAML_VERSION}"
         )
EXCLUDE_PATTERN=${EXCLUDES[0]}
for (( i=1;i<${#EXCLUDES[@]};i++ ))
do
	EXCLUDE_PATTERN="${EXCLUDE_PATTERN}|${EXCLUDES[i]}"
done
shopt -s extglob
rm -rvf "!(${EXCLUDE_PATTERN})"

echo ">>> After cleaning"
ls -la .
