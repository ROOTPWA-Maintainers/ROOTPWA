#!/bin/bash

if [[ -n "${1}" ]]
then
    svn propset svn:keywords "Rev Date Author" ${1}
else
    svn propset svn:keywords "Rev Date Author" *
fi
