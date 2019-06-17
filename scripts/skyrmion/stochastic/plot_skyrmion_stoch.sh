#!/bin/bash
gnuplot -e 'set xlabel "time [s]"; set ylabel "<mz>";p "'$1'" u 1:2 w l t "<mz>"' -persist
