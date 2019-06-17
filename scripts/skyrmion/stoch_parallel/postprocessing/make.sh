#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
module unload gcc
module load gcc/7.2.0
g++ -std=c++14 -o calc_mean_annihilationtime calc_mean_annihilationtime.cpp
