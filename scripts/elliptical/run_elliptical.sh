#!/bin/bash
# Usage .sh absolute/path/to/write/output/ <optional-GPU number>
set -e

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..
buildfile=main_elliptical.cpp
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
$magafdir/scripts/bash/check_write_dir.sh $1
if [ "$?" == 1 ];then
    exit 1
fi

# copying files
cp $magafdir/bin/magnum.af-* $1
cp $buildfile $1

# running
if [ -e $magafdir/bin/magnum.af-cuda ];then
    screen -d -m bash -c "$magafdir/bin/magnum.af-cuda $1 $2 > $1/cout.txt 2>&1"
elif [ -e $magafdir/bin/magnum.af-opencl ];then
    screen -d -m bash -c "$magafdir/bin/magnum.af-opencl $1 $2 > $1/cout.txt 2>&1"
else
    $magafdir/bin/magnum.af-cpu $1 $GPU
fi

# executing tail -f on cout.dat
if [ -e $magafdir/bin/magnum.af-cuda ];then
    tail -f $1/cout.txt
elif [ -e $magafdir/bin/magnum.af-opencl ];then
    tail -f $1/cout.txt
fi
