#!/bin/bash
A=50
B=5
outputdir=A"$A"mTB"$B"/
../../magnum.af -f -i stress_CoPt.py -o A"$A"mTB"$B" $A $B

cd $outputdir
gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "mz_over_hz.pdf";
    set xlabel "H_z [T]";
    set ylabel "m_z";
    set grid;
    plot "m.dat" u 8:5 w lp notitle
'
