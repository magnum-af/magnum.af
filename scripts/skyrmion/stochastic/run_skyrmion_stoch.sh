#!/bin/bash
# Usage .sh absolute/path/to/write/to

echo "Enter dt [s]"
read dt
echo "dt = $dt [s]"

echo "Enter T [K]"
read T
echo "T = $T [K]"

GPU=0

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../../..
buildfile=main_skyrmion_stoch.cpp
plotfile=plot*.sh

# checking if other main exists in /src
$magafdir/scripts/bash/check_main.sh

# building
cp $buildfile $magafdir/src
$magafdir/scripts/bash/build.sh $magafdir
rm $magafdir/src/$buildfile

# moving possible old main back
$magafdir/scripts/bash/check_main_remove.sh

# creating write dir
$magafdir/scripts/bash/check_write_dir.sh $1/skyrm

# copying files
cp $magafdir/bin/magnum.af-* $1
cp $buildfile $1
cp $plotfile $1

# running
if [ -e $magafdir/bin/magnum.af-opencl ];then
    screen -d -m bash -c "$magafdir/bin/magnum.af-opencl $1 $GPU $T $dt > $1/cout.txt"
else
    $magafdir/bin/magnum.af-cpu $1 $GPU $T $dt
fi

echo "To follow cout.dat run:"
echo "tail -f $1/cout.txt"

# run plot
#./plot_skyrmion_stoch.sh $1/m0.dat
