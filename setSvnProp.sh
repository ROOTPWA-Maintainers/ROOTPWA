#!/bin/bash

if [[ -n "${1}" ]]
then
    svn propset svn:keywords "Rev Date Author" ${1}
else
		for i in $(find . -type f -name '*')
		do
				svn propset svn:keywords "Rev Date Author" ${i}
		done
fi
