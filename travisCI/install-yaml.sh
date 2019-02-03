#!/bin/bash
set -exv

echo ">>> Using YAML version ${YAML_VERSION}"

cd "${TRAVIS_BUILD_DIR}"/deps/

if [ -d yaml-${YAML_VERSION} ]
then
	echo "    Existing YAML installation found in 'yaml-${YAML_VERSION}', using that."
else
	echo "    No YAML installation found, installing a fresh one."

	# clone the YAML git repository
	git clone https://github.com/jbeder/yaml-cpp.git yaml-${YAML_VERSION}

	cd yaml-${YAML_VERSION}

	# checkout one particular version we know to be working
	git checkout yaml-cpp-${YAML_VERSION}

	# compile
	mkdir build
	cd build
	cmake -DBUILD_SHARED_LIBS=ON ..
	make

	cd "${TRAVIS_BUILD_DIR}"/deps/
fi

ln -sfn yaml-${YAML_VERSION} yaml
