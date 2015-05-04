#!/bin/bash
set -exv
wget http://www.hyperrealm.com/libconfig/libconfig-1.4.9.tar.gz
tar -xzf libconfig-1.4.9.tar.gz
mkdir libconfig_install
export PARENT_DIR=$PWD
cd ./libconfig-1.4.9
./configure --prefix=$PARENT_DIR/libconfig_install
make && make install