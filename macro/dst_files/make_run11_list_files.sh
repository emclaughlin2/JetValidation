#!/bin/bash

#if [ -n $1 ] && [ $1 -le 10000 ] && [ $((65000%$1)) -eq 0 ]; then
lines=5
rm run11_dst_truth_jet_* run11_dst_global_* run11_dst_calo_cluster_*
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
split $args run11_dst_truth_jet.list "run11_dst_truth_jet_"
split $args run11_dst_global.list "run11_dst_global_"
split $args run11_dst_calo_cluster.list "run11_dst_calo_cluster_"
#else
#	echo "Please pick a number of files to have that splits 650,000 lines evenly that is less than or equal to 10,000"
#fi
