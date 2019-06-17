#!/bin/bash
gnuplot -e 'set terminal pdf enhanced;
            set output "m.pdf";
            set xlabel "t [ns]";
            set ylabel "average magnetization";
            set title "\muMAG Standard Problem 4";
            p "m.dat" u ($1*1e9):2 w l t "mx","m.dat" u ($1*1e9):3 w l t "my", "m.dat" u ($1*1e9):4 w l t "mz";
            set terminal pop;
            set output;
            p "m.dat" u ($1*1e9):2 w l t "mx","m.dat" u ($1*1e9):3 w l t "my", "m.dat" u ($1*1e9):4 w l t "mz";
            ' -persist
