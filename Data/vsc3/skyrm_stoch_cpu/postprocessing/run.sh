#!/bin/bash
# Run from parent folder 
#module unload gcc
#module load gcc/7.2.0
g++ -std=c++14 -o postprocessing/calc_mean_annihilationtime.exe postprocessing/calc_mean_annihilationtime.cpp 

echo $1
if [[ -z "$1" && (-z "$2") ]]; then
    echo "Usage: ./run.sh /path/to/inputfile /path/to/outputfolder"
    echo "e.g.:  ./run.sh \$PWD/annihilationtime.dat \$PWD"
    exit 1
fi
./postprocessing/calc_mean_annihilationtime.exe $1 $2
