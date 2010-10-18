#!/usr/bin/bash
##########################################################################
#
#    Copyright 2010
#
#    This file is part of rootpwa
#
#    rootpwa is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    rootpwa is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      sets up ROOTPWA environment needed for compilation and
#      execution of ROOTPWA programs and scripts
#
#      defines ROOTPWA environment variables
#      PWA_ENV_SET
#      ROOTPWA
#      ROOTPWA_BIN
#      ROOTPWA_LIB
#      BOOST_ROOT
#      LIBCONFIG
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# set ROOTPWA paths
export ROOTPWA=$(readlink --canonicalize $(dirname ${BASH_SOURCE}))
if [[ ! -d "${ROOTPWA}/build" ]]
then
    mkdir --parents --verbose ${ROOTPWA}/build
fi
export ROOTPWA_BIN=${ROOTPWA}/build/bin
export ROOTPWA_LIB=${ROOTPWA}/build/lib

# base directory
PWA_BASE_DIR=${ROOTPWA%/}  # remove trailing /
PWA_BASE_DIR=${PWA_BASE_DIR%/$(basename ${ROOTPWA})}  # one level above

# set Boost path
export BOOST_ROOT=${PWA_BASE_DIR}/../../boost

# set libconfig path
export LIBCONFIG=${PWA_BASE_DIR}/../../libconfig

# set paths
# add new paths to front to insure that they override older definitions
export PATH="${ROOTPWA_BIN}:${PATH}"
export LD_LIBRARY_PATH="${ROOTPWA_LIB}:${LIBCONFIG}/lib:${LD_LIBRARY_PATH}"

# report
echo ">>> info: setting ROOTPWA environment: "
echo "        ROOTPWA      = '${ROOTPWA}'"
echo "        ROOTPWA_BIN  = '${ROOTPWA_BIN}'"
echo "        ROOTPWA_LIB  = '${ROOTPWA_LIB}'"
echo "        BOOST_ROOT   = '${BOOST_ROOT}'"
echo "        LIBCONFIG    = '${LIBCONFIG}'"

# set flag that indicates that environment was set successfully
export PWA_ENV_SET="TRUE"
