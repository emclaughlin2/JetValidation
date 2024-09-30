#!/bin/bash

#if [ -n $1 ] && [ $1 -le 10000 ] && [ $((65000%$1)) -eq 0 ]; then
lines=5
rm detroit_jet10_dst_calo_waveform_* detroit_jet10_g4hits_*
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
split $args detroit_jet10_dst_calo_waveform.list "detroit_jet10_dst_calo_waveform_"
split $args detroit_jet10_g4hits.list "detroit_jet10_g4hits_"
#else
#	echo "Please pick a number of files to have that splits 650,000 lines evenly that is less than or equal to 10,000"
#fi
