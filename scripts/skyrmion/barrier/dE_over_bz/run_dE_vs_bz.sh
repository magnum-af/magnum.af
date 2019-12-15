#!/bin/bash
#usage ./script path/

bz_init=0.0
bz_step=0.005
demag=1 # 0 is nodemag, else is demag
for i in {0..10..1}
do
    bz=$(awk "BEGIN {print $bz_init+$i*$bz_step; exit}")
    echo "i=$i, bz=$bz"
    ../../../magnum.af main_bloch_zeeman_ref_Jex.cpp $1/$bz $bz $demag
done
