#!/bin/bash
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
#      collection of useful bash functions
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# converts relative into absolute pathnames
# usage: pathRel2Abs "${RELATIVE_PATH}" ABSOLUTE_PATH
function pathRel2Abs ()
{
    local _ABSOLUTE_PATH=$(readlink --canonicalize-missing ${1})
    eval "$2=${_ABSOLUTE_PATH}"
}


# prepends a path to path list stored in an environment variable making
# sure not to double entries; seperator ':' is assumed
# usage: prependPathToEnvVar VARIABLE "${PATH}"
function prependPathToEnvVar ()
{
    local _PATH=${2}
    local _VARIABLE=${1}
    local _SEPARATOR=":"
    if eval "test -z \$${_VARIABLE}"
    then
				export ${_VARIABLE}=${_PATH}
    else
				local _FOUND=$(eval echo \$${_VARIABLE} | tr ${_SEPARATOR} \\n | grep ^${_PATH}\$)
				if [[ -z "${_FOUND}" ]]
				then
            export ${_VARIABLE}=$(eval echo ${_PATH}${_SEPARATOR}\$${_VARIABLE})
				fi
    fi
}


# appends a path to path list stored in an environment variable making
# sure not to double entries; seperator ':' is assumed
# usage: appendPathToEnvVar VARIABLE "${PATH}"
function appendPathToEnvVar ()
{
    local _PATH=${2}
    local _VARIABLE=${1}
    local _SEPARATOR=":"
    if eval "test -z \$${_VARIABLE}"
    then
				export ${_VARIABLE}=${_PATH}
    else
				local _FOUND=$(eval echo \$${_VARIABLE} | tr ${_SEPARATOR} \\n | grep ^${_PATH}\$)
				if [[ -z "${_FOUND}" ]]
				then
            export ${_VARIABLE}=$(eval echo \$${_VARIABLE}${_SEPARATOR}${_PATH})
				fi
    fi
}


# waits until number of concurrently running processes is below threshold
# usage: waitForJobs "${JOB_PATTERN}" ["${MAX_NMB_JOBS}"] ["${POLL_INTERVAL}"] ["${GUARD_INTERVAL}"]
function waitForJobs ()
{
    local _JOB_PATTERN=${1}
    local _MAX_NMB_JOBS=${2:-3}
    local _POLL_INTERVAL=${3:-10s}
    local _GUARD_INTERVAL=${4:-5s}
    sleep ${_GUARD_INTERVAL}  # give job some time to show up in process list
    # infinite loop
    while true
    do
        # get number of running jobs
        declare -i _NUM_JOBS=$(pgrep -f "${_JOB_PATTERN}" | wc -l)
        # wait until number of jobs is smaller than maximum allowed
        if (( "${_NUM_JOBS}" < "${_MAX_NMB_JOBS}" ))
				then
						break
				else
						sleep ${_POLL_INTERVAL}
				fi
    done
}


# checks whether variable content is number
# usage: isNumber "${VARIABLE}"
function isNumber ()
{
    if [ ${1} -eq ${1} 2> /dev/null ]
    then
				return 0
    else
				return 1
    fi
}


# converts string(s) passed as argument(s) to lower case
# usage: newvar=$( toLower "${oldVar}" )
function toLower ()
{
    if [ -z "${1}" ]  # no argument
    then
				return
    fi
    # lower case all arguments
    echo "$@" | tr '[:upper:]' '[:lower:]'
    return
}


# converts string(s) passed as argument(s) to upper case
# usage: newvar=$( toUpper "${oldVar}" )
function toUpper ()
{
    if [ -z "${1}" ]  # no argument
    then
				return
    fi
    # upper case all arguments
    echo "$@" | tr '[:lower:]' '[:upper:]'
    return
}


# extracts value from text file
# usage: getFieldValue "${FILE_NAME}" "${VALUE_STRING}" VALUE
function getFieldValue ()
{
    local _FILE_NAME=${1}
    local _VALUE_STRING=${2}
    if [[ "${_FILE_NAME##*.}" == "gz" ]]
    then
				GREP="zgrep"
    else
				GREP="grep"
    fi
    local _VALUE_LINE=$(${GREP} "${_VALUE_STRING}" ${_FILE_NAME})
    if [[ -z "${_VALUE_LINE}" ]]
    then
        return 1
    else
				local _VALUE=${_VALUE_LINE#*${_VALUE_STRING}}
	# make sure _VALUE is number
				if [ ! ${_VALUE} -eq ${_VALUE} 2> /dev/null ]
				then
						echo "    value for '${_VALUE_STRING}' in '${_FILE_NAME}' is not an integer: ${_VALUE}"
				fi
    fi
    # deference third argument
    eval "$3=${_VALUE}"
    return 0
}
