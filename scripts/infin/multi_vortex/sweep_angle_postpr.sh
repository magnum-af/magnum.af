#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"
[ -f demag_values.dat ] && mv demag_values.dat demag_values.dat.old
for dir in $(ls -dv */); do
    echo $dir
    cat $dir/demag.dat >> demag_values.dat
    echo "" >> demag_values.dat
done

Hx_0degree=$(awk -F ', ' 'NR==1{print $4}' demag_values.dat)
echo "$Hx_0degree"

gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_demag_x_over_angle.pdf";
    set xlabel "angle [°]";
    set ylabel "Hx_{demag} [mT]";
    set grid;
    plot "demag_values.dat" u 8:($4*1e3) w lp title "simulation", "demag_values.dat" u 8:(1e3*'"$Hx_0degree"'*cos($8/180*3.141)) w lp t "Hx(phi=0) * cos(phi)"
'
gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_demag_y_over_angle.pdf";
    set xlabel "angle [°]";
    set ylabel "Hy_{demag} [mT]";
    set grid;
    plot "demag_values.dat" u 8:($5*1e3) w lp title "simulation", "demag_values.dat" u 8:(1e3*'"$Hx_0degree"'*sin($8/180*3.141)) w lp t "Hx(phi=0) * sin(p    hi)"
'
gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "h_norm_over_angle.pdf";
    set xlabel "angle [°]";
    set ylabel "|H_{demag}| [mT]";
    set grid;
    plot "demag_values.dat" u 8:(sqrt(($4*1e3)**2+($5*1e3)**2+($6*1e3)**2)) w lp notitle
'
