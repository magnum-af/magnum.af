#!/bin/bash
if [ -f bz_e_barriers.dat ]
then
    mv bz_e_barriers.dat prev_bz_e_barriers.dat
fi
for file in `ls -d */bz_over_J.dat | sort -n`; do echo $file && cat $file >> bz_e_barriers.dat; done

gnuplot plot.gnuplot

evince dE_over_bz.pdf
