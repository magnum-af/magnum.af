#!/bin/bash
gnuplot -e "
    set xlabel 't';
    p 'm.dat' u 1:2 w l, 'm.dat' u 1:3 w l, 'm.dat' u 1:4 w l, 'm.dat' u 1:(\$5*4*3.14*1e-7) w l;
" --persist

gnuplot -e "
    set xlabel 'mu0 Hx';
    set ylabel 'm_x';
    p 'm.dat' u (\$5*4*3.14*1e-7):2 w l
" --persist
