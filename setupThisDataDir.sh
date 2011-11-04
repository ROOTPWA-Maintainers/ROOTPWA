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
#      template script that sets up environment variables that are
#      used by ROOTPWA scripts to process data; please create a copy
#      and adapt it to your directory layout
#
#      relies on ROOTPWA environment variables
#      ROOTPWA_ENV_SET
#      ROOTPWA
#
#      defines ROOTPWA environment variables
#      ROOTPWA_DATA_ENV_SET  # value 'true' indicates that environment was set
#      ROOTPWA_DATA_DIR      # directory that holds mass bins
#      ROOTPWA_KEYS_DIR      # name pattern for directories with key files
#      ROOTPWA_WAVE_LIST     # path to wave list
#      ROOTPWA_FITS_DIR      # directory where fit result will be written
#      ROOTPWA_LOGS_DIR      # directory where fit log files will be written
#
#      the variables are set based on the location of _this_ file which is
#      assumed to be located in the directory that holds the mass bins
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


if [[ "${ROOTPWA_ENV_SET}" != "true" ]]
then
		echo "!!! error: ROOTPWA environment is not setup. cannot setup data environment. please source the ROOTPWA setup script first."
else

    # get path of this shell script
		THIS_FILE_PATH=${BASH_SOURCE[0]}

    # setup data directories
		export ROOTPWA_DATA_DIR=$(readlink --canonicalize $(dirname ${THIS_FILE_PATH}))


    #*************************************************************************
    # adapt to your analysis

		export ROOTPWA_KEYS_DIR="${ROOTPWA}/keyfiles/key3pi/SET?_new"
		export ROOTPWA_WAVE_LIST="${ROOTPWA}/keyfiles/key3pi/wavelist"
		export ROOTPWA_FITS_DIR="${ROOTPWA_DATA_DIR}/fits"
		export ROOTPWA_LOGS_DIR="${ROOTPWA_DATA_DIR}/logs"

    #
    #*************************************************************************


    # make sure directories for fit results and log files exist
		if [[ ! -d ${ROOTPWA_FITS_DIR} ]]
		then
				mkdir --parents --verbose ${ROOTPWA_FITS_DIR}
		fi
		if [[ ! -d ${ROOTPWA_LOGS_DIR} ]]
		then
				mkdir --parents --verbose ${ROOTPWA_LOGS_DIR}
		fi

    # check paths
		PATHS_TO_CHECK=( ${ROOTPWA_DATA_DIR} ${ROOTPWA_KEYS_DIR} ${ROOTPWA_FITS_DIR} ${ROOTPWA_LOGS_DIR} )
		for (( IDX=0; IDX<"${#PATHS_TO_CHECK[@]}"; IDX+=1 ))
		do
				if [[ ! -d ${PATHS_TO_CHECK[IDX]} ]]
				then
						echo "??? warning: directory '${PATHS_TO_CHECK[IDX]}' does not exist. ROOTPWA scripts may not work as expected."
				fi
		done
		PATHS_TO_CHECK=( ${ROOTPWA_WAVE_LIST} )
		for (( IDX=0; IDX<"${#PATHS_TO_CHECK[@]}"; IDX+=1 ))
		do
				if [[ ! -f ${PATHS_TO_CHECK[IDX]} ]]
				then
						echo "??? warning: file '${PATHS_TO_CHECK[IDX]}' does not exist. ROOTPWA scripts may not work as expected."
				fi
		done

    # set flag that indicates that environment was set successfully
		export ROOTPWA_DATA_ENV_SET="true"

		echo ">>> info: setting ROOTPWA data environment: "
		echo "        ROOTPWA_DATA_ENV_SET = '${ROOTPWA_DATA_ENV_SET}'"
		echo "        ROOTPWA_DATA_DIR     = '${ROOTPWA_DATA_DIR}'"
		echo "        ROOTPWA_KEYS_DIR     = '${ROOTPWA_KEYS_DIR}'"
		echo "        ROOTPWA_WAVE_LIST    = '${ROOTPWA_WAVE_LIST}'"
		echo "        ROOTPWA_FITS_DIR     = '${ROOTPWA_FITS_DIR}'"
		echo "        ROOTPWA_LOGS_DIR     = '${ROOTPWA_LOGS_DIR}'"
fi