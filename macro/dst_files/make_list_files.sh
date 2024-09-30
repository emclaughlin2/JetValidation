#!/bin/bash

#if [ -n $1 ] && [ $1 -le 10000 ] && [ $((65000%$1)) -eq 0 ]; then
lines=5
rm dst_truth_jet_* dst_calo_waveform_* dst_global_* dst_calo_cluster_* dst_calo_nozero_* dst_truthinfo_*
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
split $args dst_truth_jet.list "dst_truth_jet_"
split $args dst_calo_waveform.list "dst_calo_waveform_"
split $args dst_global.list "dst_global_"
split $args dst_calo_cluster.list "dst_calo_cluster_"
split $args dst_calo_nozero.list "dst_calo_nozero_"
split $args dst_truthinfo.list "dst_truthinfo_"
#else
#	echo "Please pick a number of files to have that splits 650,000 lines evenly that is less than or equal to 10,000"
#fi
