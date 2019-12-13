#!/bin/bash
gnuplot -e '
    set terminal pdf;
    set output "rigid_moment.pdf";
    set xlabel "time [s]";
    p "m.dat" u 1:2 w l t "mx", "m.dat" u 1:3 w l t "my"
'
