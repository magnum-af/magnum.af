#!/bin/bash

echo "Usage: ./run*.sh {$}PWD/path/to/write/oputput"
echo "Building main for vsc3 \n NOTE: remove opencl in src/CMakeLists.txt"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
magafdir=../../..

ismainmoved=1
# aborting if src/main* already exists
if [ -e $magafdir/src/main*.cpp ]
then
    echo "Temoraryly moving current main in /src"
    ismainmoved=0
    mv $magafdir/src/main*.cpp $magafdir
	ls $magafdir
else
    echo "src/ is clean, building..."
fi

# determining write directory
if [ -n "$1" ];then
    if [ -e "$1" ];then
        echo "Error: Write Directory exists!" $1
	exit 1
    else
        mkdir --parents $1/data/vti
        echo "writing in new directory " $1
    fi
fi

cp main_skyrmion_stoch_parallel.cpp  $magafdir/src
# Compiling 
$magafdir/scripts/vsc3/compile.sh

# cleaning up main(s)
rm $magafdir/src/main_skyrmion_stoch_parallel.cpp 
if [ $ismainmoved ]
then
    mv $magafdir/main*.cpp $magafdir/src/
fi

cp main_skyrmion_stoch_parallel.cpp  $1
cp $magafdir/bin/* $1

echo "Enter dt [s]"
read dt
echo $dt

echo "Enter T [k]"
read T
echo $T

echo "Enter Numer of runs"
read runs
echo "runs = $runs"

echo "dt=$dt \n  T=$T \n runs=$runs" >> $1/inputvars.txt
for ((i = 1; i <= $runs; i++)); do
    echo "$1/magnum.af-cpu $1 $DIR/data $dt $T $i 0" >> $1/parallel_commands.txt
done

# running slurm script
cd $1
sbatch $DIR/vsc-parallel.sh
#To test run: 
#$DIR/vsc-parallel.sh
