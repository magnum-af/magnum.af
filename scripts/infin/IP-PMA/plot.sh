#!/bin/bash
gnuplot -e 'set terminal pdf enhanced;
            set output "err_over_phi.pdf";
            set xlabel "phi [rad]";
            set ylabel "err";
            set title "Angle error of demagfield (ideal sensor)";
            p "h.dat" u 1:2 w l t "h_err","h.dat" u 1:4 w l t "m_err" ;
            set terminal pop;
            set output;
            p "h.dat" u 1:2 w l t "h_err","h.dat" u 1:4 w l t "m_err" ;
            ' -persist
