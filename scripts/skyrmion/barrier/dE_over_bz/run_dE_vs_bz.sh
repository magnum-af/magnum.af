#!/bin/bash

bz_init=0.0
bz_step=0.005
for i in {0..10..1}
do
    bz=$(awk "BEGIN {print $bz_init+$i*$bz_step; exit}")
    echo "i=$i, bz=$bz"
    ../../../magnum.af -n main_bloch_zeeman_ref_Jex_demag.cpp $1/$bz 3 main_bloch_zeeman_ref_Jex_demag.cpp $bz
done
