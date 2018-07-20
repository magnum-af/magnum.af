#!/bin/bash
# Note: to plot temperature dependent, edit mean_annihilationtimes.dat manually
module unload gcc
module load gcc/4.9 gnuplot/5.0.5

gnuplot -e 'set terminal pdf; set output "'$1'/mean_annihilationtime.pdf";set xrange[280:520];set xlabel "T [K]"; set ylabel "t [s]";plot "'$1'/mean_annihilationtimes.dat" u 5:1:2 with yerrorbars t "mean annihilation time"' -persist
gnuplot -e 'set xrange[280:520];set xlabel "T [K]"; set ylabel "t [s]";plot "'$1'/mean_annihilationtimes.dat" u 5:1:2 with yerrorbars t "mean annihilation time"' -persist
