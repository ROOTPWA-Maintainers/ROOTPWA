#!/bin/bash
set -exv
curl -L -O https://cmake.org/files/v3.2/cmake-3.2.2-Linux-x86_64.tar.gz
tar -xzf cmake-3.2.2-Linux-x86_64.tar.gz
mv ./cmake-3.2.2-Linux-x86_64 ./cmake
