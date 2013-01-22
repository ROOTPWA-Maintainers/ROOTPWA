#!/bin/bash

if [ -z $ROOTPWA ]; then
	echo '$ROOTPWA not set, set it correctly for this script to work.'
	exit 1
fi

ln -isv $ROOTPWA/gitUtils/clientHook.py $ROOTPWA/.git/hooks/commit-msg

