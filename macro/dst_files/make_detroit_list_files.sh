#!/bin/bash

#if [ -n $1 ] && [ $1 -le 10000 ] && [ $((65000%$1)) -eq 0 ]; then
lines=5
rm detroit_dst_truth_jet_* detroit_dst_calo_waveform_* detroit_dst_global_* detroit_dst_calo_cluster_* detroit_dst_calo_nozero_* detroit_dst_truthinfo_*
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
split $args detroit_dst_truth_jet.list "detroit_dst_truth_jet_"
split $args detroit_dst_calo_waveform.list "detroit_dst_calo_waveform_"
split $args detroit_dst_global.list "detroit_dst_global_"
split $args detroit_dst_calo_cluster.list "detroit_dst_calo_cluster_"
split $args detroit_dst_calo_nozero.list "detroit_dst_calo_nozero_"
split $args detroit_dst_truthinfo.list "detroit_dst_truthinfo_"
#else
#	echo "Please pick a number of files to have that splits 650,000 lines evenly that is less than or equal to 10,000"
#fi
