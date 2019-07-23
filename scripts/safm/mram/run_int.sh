#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

xyz_dir=2 # z dir
int_over_min="int" # integrate, instead of minimizing
hzee_max="0.12" # maximum H_external for hysteresis loop
writedir="$HOME"/data_magnum.af/safm/mram
../../magnum.af -f mram.cpp "$writedir" "$xyz_dir" "$int_over_min" "$hzee_max"
# -p plot.gpi

cd "$writedir"
gnuplot -e '
    set terminal pdfcairo enhanced;
    set output "mram_hys.pdf";
    set xlabel "H_ext [T]";
    set ylabel "<m_z>";
    set grid;
    set title "MRAM cell hard axis loop";
    p "m_int.dat" u 7:2 w lp t "<mx>", "" u 7:3 w lp t "<my>", "" u 7:4 w lp t "<mz>"
'
    #minimizer#p "m.dat" u 7:2 w lp t "<mx>", "" u 7:3 w lp t "<my>", "" u 7:4 w lp t "<mz>"
