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
#      runs a program with various valgrind tools
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


echo ">>> info: called ${0} ${*}"


#VALGRIND_TOOL="memcheck"
VALGRIND_TOOL="callgrind"

PROGRAM="../../build/bin/calcAmplitudes"
PROGRAM_OPT="-n 200000 -k 1-4++1+f21270_32_pi-.key -p ./particleDataTable.txt -o 1-4++1+f21270_32_pi-.amp /dev/shm/allBins.ps.noTarget.root"
#PROGRAM="../debugBuild/bin/testLibppFunctions"
#PROGRAM_OPT=""

LOG_FILE="./"$(basename ${PROGRAM})"_${VALGRIND_TOOL}.log"
OUT_FILE="./"$(basename ${PROGRAM})"_${VALGRIND_TOOL}.out"


echo ">>> info: ${0} started on $(date)"
echo ">>> info: running valgrind tool ${VALGRIND_TOOL} on '${PROGRAM}'"
echo

# prepare command line
case ${VALGRIND_TOOL} in
    "memcheck")  # run memcheck tool
	VALGRIND_OPT="--leak-check=full --suppressions=${ROOTSYS}/etc/valgrind-root.supp --show-reachable=yes --track-origins=yes --freelist-vol=100000000 --verbose"
	;;
    "callgrind")  # run callgrind tool
	#!!! callgrind segfaults when demangling Boost symbols
	#!!! workaround by disabling demangling and running c++filt as an afterburner
	#!!! see http://bugs.kde.org/show_bug.cgi?id=197988
	VALGRIND_OPT="--tool=callgrind --simulate-cache=yes --cacheuse=yes --demangle=no --callgrind-out-file=${OUT_FILE} --verbose"
	;;
    *)
	echo "!!! error: unknown valgrind tool ${VALGRIND_TOOL}. exiting."
	exit 1
	;;
esac

# do it
CMD="valgrind ${VALGRIND_OPT} ${PROGRAM} ${PROGRAM_OPT} &> ${LOG_FILE}"
echo "${CMD}"
time eval ${CMD}
echo

# demangle symbols in callgrind output
if [[ ${VALGRIND_TOOL} = "callgrind" && -f "${OUT_FILE}" ]]
then
    c++filt < ${OUT_FILE} > ${OUT_FILE}.demangled
fi

echo ">>> info: ${0} finished on $(date)"
exit 0
