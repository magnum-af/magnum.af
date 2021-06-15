#!/bin/bash

odir="output_true_pbc_data"
mkdir -p "$odir"
for N in 40 64 100 130 208
do
    echo "$odir"/run_N"$N".out
    ./true_pbc.py -o "$odir" -d 0 "$N" |& tee "$odir"/run_N"$N".out
done 
