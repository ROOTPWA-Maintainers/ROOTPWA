#!/bin/bash

for i in *.key ;
  do cp $i $i.asc;
     sed 's/ascii/binary/g' $i.asc > $i;
done ;