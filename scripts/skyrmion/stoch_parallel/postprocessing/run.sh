#!/bin/bash
# usage: ./run /path/to/simulations/directories

cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
./make.sh

echo $1
if [ -z "$1" ]; then
    echo "Usage: ./run.sh /path/to/inputdirectories"
    echo "e.g.:  ./run.sh \$PWD/annihilationtime.dat"
    exit 1
fi

for dir in $1/*
do 
    ./calc_mean_annihilationtime $dir/anihilationtime.dat $dir
done

if [ -f $1/mean_annihilationtimes.dat ]; then
    rm $1/mean_annihilationtimes.dat
fi
for dir in $1/*/
do 
    echo "#"$dir >>  $1/mean_annihilationtimes.dat
    cat $dir/mean_annihilationtime.dat >> $1/mean_annihilationtimes.dat
done
