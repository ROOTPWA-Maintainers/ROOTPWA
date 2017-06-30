#!/bin/bash
set -exv

if [ ${TRAVIS_RPWA_RUN_FORMAT:-0} -eq 1 ] ; then
	travisCI/run_format.sh
fi
if [ ${TRAVIS_RPWA_RUN_BUILD:-0} -eq 1 ] ; then
	travisCI/run_build.sh
fi
