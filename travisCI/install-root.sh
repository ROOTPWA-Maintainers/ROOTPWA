#!/bin/bash
set -exv

ROOT_VERSION=6.06.06

cd ${TRAVIS_BUILD_DIR}/deps/

if [ -d root-${ROOT_VERSION} ] ; then
	echo "ROOT installation found in 'root-${ROOT_VERSION}', using that."
else
	echo "No ROOT installation found, installing a fresh one."

	wget https://root.cern.ch/download/root_v${ROOT_VERSION}.Linux-ubuntu14-x86_64-gcc4.8.tar.gz -O root.tar.gz

        mkdir root-${ROOT_VERSION}
	tar -xzf root.tar.gz -C root-${ROOT_VERSION} --strip-components=1
	rm -rf root.tar.gz
fi
ln -sfn root-${ROOT_VERSION} root
