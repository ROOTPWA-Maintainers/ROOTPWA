#!/bin/bash
set -exv

echo ">>> Using BAT version ${BAT_VERSION}"

cd "${TRAVIS_BUILD_DIR}"/deps/

if [ -d bat-${BAT_VERSION} ]
then
	echo "    Existing BAT installation found in 'bat-${BAT_VERSION}', using that."
else
	echo "    No BAT installation found, installing a fresh one."

	# clone the BAT git repository
	git clone https://github.com/bat/bat.git bat-${BAT_VERSION}

	cd bat-${BAT_VERSION}

	# checkout one particular version we know to be working
	git checkout v${BAT_VERSION}

	# create configure and build scripts
	./autogen.sh

	# configure to be installed in the current (source) directory
	./configure --prefix=$(pwd -P) --enable-parallel

	# compile and install
	make --jobs=${NMB_JOBS}
	make install

	cd "${TRAVIS_BUILD_DIR}"/deps/
fi

rm -vf bat
ln -sfn bat-${BAT_VERSION} bat
