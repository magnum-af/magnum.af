#!/bin/bash
../magnum.af -fi sp4.cpp -o sp4
../magnum.af -fi sp4_matmul.cpp -o sp4_matmul
../magnum.af -fi sp4_matmul_nobc.cpp -o sp4_matmul_nobc -p plot_compare.sh 
