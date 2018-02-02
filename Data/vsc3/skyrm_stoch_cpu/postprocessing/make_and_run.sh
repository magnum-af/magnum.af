#!/bin/bash
module unload gcc
module load gcc/7.2.0
#gcc -v
g++ -std=c++14 -o calc_mean_annihilationtime.exe calc_mean_annihilationtime.cpp 
./calc_mean_annihilationtime.exe $PWD/anihilationtime.dat $PWD

