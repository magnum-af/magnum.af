#!/bin/bash
if [ ! -n "$1" ]; then
    echo "Usage: ./run*.sh /absolute/path/to/write/to"
    exit 1
fi

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
magafdir=../..
buildfile=main_sp4_newllg.cpp
plotfile=plot_sp4.sh

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
    $magafdir/bin/magnum.af-opencl $1
elif [ -e $magafdir/bin/magnum.af-cuda ];then
    $magafdir/bin/magnum.af-cuda $1
else
    $magafdir/bin/magnum.af-cpu $1
fi

# plotting
./$plotfile $1/m.dat
