#!/bin/bash
set -exv

echo ">>> Using BAT revision ${BAT_REVISION}"

cd "${TRAVIS_BUILD_DIR}"/deps/

if [ -d bat-${BAT_REVISION} ]
then
	echo "    Existing BAT installation found in 'bat-${BAT_REVISION}', using that."
else
	echo "    No BAT installation found, installing a fresh one."

	# clone the BAT git repository
	git clone https://github.com/bat/bat.git bat-${BAT_REVISION}

	cd bat-${BAT_REVISION}

	# checkout one particular version we know to be working
	git checkout ${BAT_REVISION}

	# create configure and build scripts
	./autogen.sh

	# configure to be installed in the current (source) directory
	./configure --prefix=$PWD --enable-parallel

	# compile and install
	make
	make install

	cd "${TRAVIS_BUILD_DIR}"/deps/
fi

ln -sfn bat-${BAT_REVISION} bat
