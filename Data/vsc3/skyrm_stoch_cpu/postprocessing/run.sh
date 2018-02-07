#!/bin/bash
# Run from parent folder 
#module unload gcc
#module load gcc/7.2.0
g++ -std=c++14 -o postprocessing/calc_mean_annihilationtime.exe postprocessing/calc_mean_annihilationtime.cpp 

echo $1
if [[ ! -e "$1" && (! -d "$2") ]]; then
    echo "Usage: ./run.sh /path/to/inputfile /path/to/outputfolder"
    echo "e.g.:  ./run.sh \$PWD/annihilationtime.dat \$PWD"
else
    ./postprocessing/calc_mean_annihilationtime.exe $1 $2
fi
