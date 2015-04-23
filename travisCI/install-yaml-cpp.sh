#!/bin/bash
set -exv
wget https://github.com/jbeder/yaml-cpp/archive/release-0.5.2.tar.gz
tar -xzf release-0.5.2.tar.gz
export PARENT_DIR=${PWD}
mkdir -p ${PARENT_DIR}/yaml-cpp-build
cd ${PARENT_DIR}/yaml-cpp-build
cmake ${PARENT_DIR}/yaml-cpp-release-0.5.2 -DCMAKE_INSTALL_PREFIX=${PARENT_DIR}/yaml-cpp-install -DBUILD_SHARED_LIBS=ON
make && make install
