#!/bin/bash
set -exv

echo ">>> System information:"
uname --all
lscpu
lsblk
df --human-readable
cat /proc/meminfo
free --mega
