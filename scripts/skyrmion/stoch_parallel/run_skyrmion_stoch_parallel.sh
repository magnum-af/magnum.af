#!/bin/bash

# run with eg:
# ~/git/magnum.af/scripts/skyrmion/stoch_parallel/run_skyrmion_stoch_parallel.sh $PWD/skyrm_stoch

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
magafdir=../../..

# checking if src/main* already exists
ismainmoved=1
if [ -e $magafdir/src/main*.cpp ]
then
    echo "Temoraryly moving current main in /src to .."
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
        mkdir --parents $1/data
        echo "writing in new directory " $1
    fi
fi

# building
cp main_skyrmion_stoch_parallel.cpp  $magafdir/src
$magafdir/scripts/build.sh $magafdir

# cleaning up main(s)
rm $magafdir/src/main_skyrmion_stoch_parallel.cpp 
if [ $ismainmoved ]
then
    mv $magafdir/main*.cpp $magafdir/src/
fi

# copying into $1
cp main_skyrmion_stoch_parallel.cpp  $1
cp $magafdir/bin/* $1

# enter parameters
echo "Enter dt [s]"
read dt
echo "dt = $dt [s]"

echo "Enter T [K]"
read T
echo "T = $T [K]"

echo "Enter Numer of runs"
read runs
echo "runs = $runs"

# writing run-file for parallel
for ((i = 1; i <= $runs; i++)); do
    echo $i
    echo $dt
    echo $T
    echo "$1/magnum.af-cpu $1 $DIR/data $dt $T $i 0" >> $1/parallel_commands.txt
done

# running
cd $1
echo "Starting Simulation in screen"
screen -d -m bash -c "parallel < parallel_commands.txt --no-notice"
