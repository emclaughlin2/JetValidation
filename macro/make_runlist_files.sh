#!/bin/bash

#if [ -n $1 ] && [ $1 -le 10000 ] && [ $((65000%$1)) -eq 0 ]; then
lines=5000
rm runlist_* 
args="-l $lines --numeric-suffixes=0 --suffix-length=1 --additional-suffix=.txt"
split $args runlist.txt "runlist_"
#else
#	echo "Please pick a number of files to have that splits 650,000 lines evenly that is less than or equal to 10,000"
#fi
