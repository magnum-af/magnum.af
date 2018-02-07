#!/bin/bash
#echo $PWD
for dir in ./*/
do
    ./postprocessing/run.sh ${dir}/data/anihilationtime.dat $PWD
done

#unload module gcc
#module load gcc/5.3 gnuplot/5.0.5
gnuplot -e 'set xrange[0:1e-13]; p "mean_annihilationtime.dat" u 1:2:3 with errorbars' -persist
gnuplot -e 'set terminal pdfcairo; set output "mean_annihilationtime.pdf"; set xrange[0:1e-13]; p "mean_annihilationtime.dat" u 1:2:3 with errorbars'
