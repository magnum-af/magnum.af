#!/bin/bash
#$1 path to write to (ending with "/")

for i in {10..20..10}; do
    echo $i
    ../../../magnum.af ../multi_vortex_demag.py $1$i $i
done

cd $1
for dir in $(ls -dv */); do
    echo $dir
    cat $dir/demag.dat >> demag_values.dat
    echo "" >> demag_values.dat
done

gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_demag_x_over_nx.pdf";
    set xlabel "nx";
    set ylabel "Hx_{demag} [mT]";
    set title "x-component of H_{demag} at (0,0)";
    plot "demag_values.dat" u 1:($3*1e3) w lp notitle
'
