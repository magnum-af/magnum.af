#!/bin/bash

magafdir=../..
# aborting if src/main* already exists
if [ -e $magafdir/src/main*.cpp ]
then
    echo "To run, please remove main in magnun.af/src! Aborting..."
    exit 1
else
    echo "src/ is clean, building..."
fi

# determining write directory
if [ -n "$1" ];then
    if [ -e "$1" ];then
        echo "writing in existing directory " $1
    else
        mkdir --parents $1
        echo "writing in new directory " $1
    fi
else
    if [ ! -e "$magafdir/Data/Testing" ];then
        echo " creating default directory $magafdir/Data/Testing"
        mkdir --parents $magafdir/Data/Testing
    else
        echo "writing in existing default directory $magafdir/Data/Testing"
    fi
fi

# building
cp main_sp4.cpp $magafdir/src
$magafdir/scripts/build.sh
rm $magafdir/src/main_sp4.cpp

# running
if [ -e $magafdir/bin/magnum.af-opencl ];then
    $magafdir/bin/magnum.af-opencl $1
elif [ -e $magafdir/bin/magnum.af-cuda ];then
    $magafdir/bin/magnum.af-cuda $1
else
    $magafdir/bin/magnum.af-cpu $1
fi

# plotting
if [ -n "$1" ]; then
    ./plot_sp4.sh $1/m.dat
else 
    ./plot_sp4.sh $magafdir/Data/Testing/m.dat
fi
