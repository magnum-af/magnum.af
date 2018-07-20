#!/bin/bash
#if [ -z "$1" ]; then
if [ $# -eq 0 ]
  then
    echo "Usage: ./plot.sh /path/to/inputfiles"
    echo "e.g.:  ./plot.sh \$PWD/..."
    exit 1
fi

echo "enter timestep"
read timestep
echo "enter maxID"
read maxID
gnuplot -e 'plot for [i=0:'$maxID'] "'$1'/m".i.".dat"     using ($0*'$timestep'):1 title "mz".i' -persist
gnuplot -e 'plot for [i=0:'$maxID'] "'$1'/m_avg".i.".dat" using ($0*'$timestep'):1 title "mz_avg".i' -persist
