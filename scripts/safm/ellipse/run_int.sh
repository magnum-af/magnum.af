#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

writedir="$HOME"/data_magnum.af/safm/ellipse/
xyz_dir=2 # z dir
int_over_min="int" # integrate, instead of minimizing
hzee_max="2.5" # maximum H_external for hysteresis loop
integr_time="100e-10"
../../magnum.af -f elipse_vortex_prop.cpp "$writedir" "$xyz_dir" "$int_over_min" "$hzee_max" "$integr_time"
# -p plot.gpi

cd "$writedir"
gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "ellipse_hys.pdf";
    set xlabel "H_{ext} [T]";
    set ylabel "<m_z>";
    set grid;
    set title "Ellipse hard axis loop";
    p "m.dat" u 7:4 w l t "<mz>";
'
    #p "m.dat" u 7:2 w lp t "<mx>", "" u 7:3 w lp t "<my>", "" u 7:4 w lp t "<mz>"
    #minimizer#p "m.dat" u 7:2 w lp t "<mx>", "" u 7:3 w lp t "<my>", "" u 7:4 w lp t "<mz>"
