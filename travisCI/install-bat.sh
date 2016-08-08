#!/bin/bash
set -exv

# clone the BAT git repository
git clone https://github.com/bat/bat.git

cd bat

# checkout one particular version we know to be working
# this version has been checked 08.08.2016
git checkout bede80ec026b381e12b5436f7cc05d12c1ec07c7

# create configure and build scripts
./autogen.sh

# configure to be installed in the current (source) directory
./configure --prefix=$PWD --enable-parallel

# compile and install
make
make install
