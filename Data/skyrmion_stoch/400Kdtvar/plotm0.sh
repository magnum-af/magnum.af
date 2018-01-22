#!/bin/bash
for dir in $PWD/*/
do
    dir=${dir%*/}
    cd ${dir}
    echo ${dir##*/}
    gnuplot -p -e 'set format x "%.3e"; p "m0.dat" u 1:4 t "mz"'
    gnuplot -p -e 'set format x "%.3e"; p "m0.dat" u 1:5 t "avg"'
    #tail -n 10000 m0.dat | gnuplot -p -e '...'
done
