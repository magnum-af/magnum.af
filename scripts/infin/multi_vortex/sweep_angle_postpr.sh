#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"
[ -f demag_values.dat ] && mv demag_values.dat demag_values.dat.old
for dir in $(ls -dv */); do
    echo $dir
    cat $dir/demag.dat >> demag_values.dat
    echo "" >> demag_values.dat
done

gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_demag_x_over_angle.pdf";
    set xlabel "angle [°]";
    set ylabel "Hx_{demag} [mT]";
    set grid;
    plot "demag_values.dat" u 8:($4*1e3) w lp notitle
'
