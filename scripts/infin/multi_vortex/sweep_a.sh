#!/bin/bash
#$1 path to write to (ending with "/")

for a in {500..2000..100}; do # a is the spacing in nm between disks
    echo $a
    ../../magnum.af multi_vortex_demag.py $1$a 100 $a
done

cd $1
for dir in $(ls -dv */); do
    echo $dir
    cat $dir/demag.dat >> demag_values.dat
    echo "" >> demag_values.dat
done

gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_demag_x_over_a.pdf";
    set xlabel "a [nm]";
    set ylabel "Hx_{demag} [mT]";
    set grid;
    set title "x-component of H_{demag} at (0,0)";
    plot "demag_values.dat" u 3:($4*1e3) w lp notitle
'
