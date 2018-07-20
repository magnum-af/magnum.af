#!/bin/bash
#if [ -z "$1" ]; then
if [ $# -eq 0 ]
  then
    echo "Usage: ./plot.sh /path/to/m*.dat"
    exit 1
fi

module unload gcc
module load gcc/4.9 gnuplot/5.0.5
echo "enter timestep"
read timestep
gnuplot -e 'set terminal pdf; set output "'$1'.pdf";set xlabel "t [s]"; set ylabel "<m_z>";plot "'$1'" using ($0*'$timestep'):1 title "<m_z>"'
