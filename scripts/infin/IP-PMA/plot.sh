#!/bin/bash
gnuplot -e 'set terminal pdf enhanced;
            set output "err_over_phi.pdf";
            set xlabel "phi [degree]";
            set ylabel "err";
            set title "Angle error of demagfield (ideal sensor)";
            p "h.dat" u 0:7 w l t "err","h.dat" u 0:8 w l t "err" ;
            set terminal pop;
            set output;
            p "h.dat" u 0:7 w l t "err","h.dat" u 0:8 w l t "err" ;
            ' -persist
