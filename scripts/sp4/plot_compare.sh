#!/bin/bash
gnuplot -e 'set terminal pdf enhanced; 
            set output "m.pdf"; 
            set xlabel "t [ns]";
            set ylabel "average magnetization";
            set title "\muMAG Standard Problem 4";
            p "m.dat" u ($1*1e9):3 w l t "sparcenobc my", "../sp4_matmul/m.dat" u ($1*1e9):3 w l t "matmul my", "../sp4/m.dat" u ($1*1e9):3 w l t "conv";
            set terminal pop;
            set output;
            replot;
            ' -persist
