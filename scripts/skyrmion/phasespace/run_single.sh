#!/bin/bash
# Usage:
# run.sh main.cpp absolute/path/to/write/output/ GPU-number dmi Ku1
set -e

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../../..
buildfile=$1
plotfile=plot*

# creating write dir
$magafdir/scripts/bash/check_write_dir.sh $2
if [ "$?" -ne 0 ];then
    exit 1
fi

# checking if other main exists in /src
$magafdir/scripts/bash/check_main.sh

# building
cp $buildfile $magafdir/src
$magafdir/scripts/bash/build.sh $magafdir
rm $magafdir/src/$buildfile

# moving possible old main back
$magafdir/scripts/bash/check_main_remove.sh

# copying files
cp $magafdir/bin/magnum.af-* $2
cp $buildfile $2
cp $plotfile $2

# running
if [ -e $magafdir/bin/magnum.af-cuda ];then
    screen -d -S GPU$3 -m bash -c "$magafdir/bin/magnum.af-cuda $2 $3 $4 $5 > $2/cout.txt 2>&1"
elif [ -e $magafdir/bin/magnum.af-opencl ];then
    screen -d -S GPU$3 -m bash -c "$magafdir/bin/magnum.af-opencl $2 $3 $4 $5 > $2/cout.txt 2>&1"
else
    $magafdir/bin/magnum.af-cpu $2 $3 $4 $5
fi

echo "To follow cout.dat run:"
echo "tail -f $2/cout.txt"
