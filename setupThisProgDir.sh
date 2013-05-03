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
#
# Description:
#      template script that sets up ROOTPWA environment needed for
#      compilation and execution of ROOTPWA programs and scripts;
#      please create a copy and adapt it to your directory layout
#
#      defines the following environment variables:
#      ROOTPWA_ENV_SET  # value 'true' indicates that environment was set
#      ROOTPWA          # ROOTPWA root directory
#      ROOTPWA_BIN      # directory that holds ROOTPWA binaries
#      ROOTPWA_LIB      # directory that holds ROOTPWA libraries
#      BOOST_ROOT       # path of Boost installation
#      LIBCONFIG        # path of libconfig installation
#
#      the variables are set based on the path of _this_ file which is
#      assumed to be located in the root directory of a working
#      ROOTPWA installation
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# get path of this shell script
THIS_FILE_PATH=${BASH_SOURCE[0]}

# set ROOTPWA paths
export ROOTPWA=$(readlink --canonicalize $(dirname ${THIS_FILE_PATH}))
if [[ ! -d "${ROOTPWA}/build" ]]
then
    mkdir --parents --verbose ${ROOTPWA}/build
fi
export ROOTPWA_BIN=${ROOTPWA}/build/bin
export ROOTPWA_LIB=${ROOTPWA}/build/lib


#*************************************************************************
# adapt to your directory layout

# base directory for auxiliary libraries
# in this example it is the directory one above the ROOTPWA directory
AUX_LIB_DIR=${ROOTPWA%/}  # remove trailing /
AUX_LIB_DIR=${AUX_LIB_DIR%/$(basename ${ROOTPWA})}  # one level above

# set Boost path
export BOOST_ROOT=${AUX_LIB_DIR}/boost

# set libconfig path
export LIBCONFIG=${AUX_LIB_DIR}/libconfig

#
#*************************************************************************


# set paths
# add new paths to front to insure that they override older definitions
export PATH="${ROOTPWA_BIN}:${PATH}"
export LD_LIBRARY_PATH="${ROOTPWA_LIB}:${LIBCONFIG}/lib:${LD_LIBRARY_PATH}"

# set flag that indicates that environment was set successfully
export ROOTPWA_ENV_SET="true"

# check paths
PATHS_TO_CHECK=( ${ROOTPWA} ${BOOST_ROOT} ${LIBCONFIG} ${LIBCONFIG}/lib )
for (( IDX=0; IDX<"${#PATHS_TO_CHECK[@]}"; IDX+=1 ))
do
		if [[ ! -d ${PATHS_TO_CHECK[IDX]} ]]
		then
				echo "??? warning: directory '${PATHS_TO_CHECK[IDX]}' does not exist. ROOTPWA may not work as expected."
		fi
done


# report
echo ">>> info: setting ROOTPWA environment variables: "
echo "        ROOTPWA_ENV_SET = '${ROOTPWA_ENV_SET}'"
echo "        ROOTPWA         = '${ROOTPWA}'"
echo "        ROOTPWA_BIN     = '${ROOTPWA_BIN}'"
echo "        ROOTPWA_LIB     = '${ROOTPWA_LIB}'"
echo "        BOOST_ROOT      = '${BOOST_ROOT}'"
echo "        LIBCONFIG       = '${LIBCONFIG}'"
