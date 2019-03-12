#!/bin/bash
set -exv

cd "${TRAVIS_BUILD_DIR}"

BASE_COMMIT=$(git rev-parse ${TRAVIS_BRANCH})
echo ">>> Running clang-format against branch ${TRAVIS_BRANCH} with hash ${BASE_COMMIT}"
RESULT_OUTPUT="$(git-clang-format --commit ${BASE_COMMIT} --diff)"

if [ "${RESULT_OUTPUT}" == "no modified files to format" ] || [ "${RESULT_OUTPUT}" == "clang-format did not modify any files" ]
then
	echo "*** Success: clang-format passed."
else
	echo "??? Warning: clang-format failed with the following output:"
	echo "${RESULT_OUTPUT}"
fi
