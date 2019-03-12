#!/bin/bash
gnuplot -e '
set terminal pdfcairo enhanced;
set output "m_over_t.pdf";
set xlabel "t [s]";
set ylabel "<mx>";
p "m.dat" u 1:2 w l t "<mx>", "" u 1:3 w l t "<my>", "" u 1:4 w l t "<mz>"
'
gnuplot -e '
set terminal pdfcairo enhanced;
set output "mx_over_t.pdf";
set xlabel "t [s]";
set ylabel "<mx>";
p "m.dat" u 1:2 w l t "<mx>"
'
gnuplot -e '
set terminal pdfcairo enhanced;
set output "my_over_t.pdf";
set xlabel "t [s]";
set ylabel "<my>";
p "m.dat" u 1:3 w l t "<my>"
'
gnuplot -e '
set terminal pdfcairo enhanced;
set output "mz_over_t.pdf";
set xlabel "t [s]";
set ylabel "<mz>";
p "m.dat" u 1:4 w l t "<mz>"
'
gnuplot -e '
set xlabel "t [s]";
set ylabel "<mx>";
p "m.dat" u 1:2 w lp t "<mx>", "" u 1:3 w lp t "<my>", "" u 1:4 w lp t "<mz>"
' --persist
