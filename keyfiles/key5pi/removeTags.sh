#!/bin/bash

for i in *.key;
do
A="$i";
B=${A//\(/};
C=${B//\)/};
D=${C//\[/_};
E=${D//\]/_};
F=${E//\>/=};
cp $i "STRIPED/$F";
done