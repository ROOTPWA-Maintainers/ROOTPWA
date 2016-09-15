#!/bin/bash
set -exv
wget https://cmake.org/files/v3.6/cmake-3.6.1-Linux-x86_64.tar.gz
tar -xzf cmake-3.6.1-Linux-x86_64.tar.gz
mv ./cmake-3.6.1-Linux-x86_64 ./cmake
