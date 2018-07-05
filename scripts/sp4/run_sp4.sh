#!/bin/bash

if [ -e ../../src/main*.cpp ]
then
    echo "To run, please remove main in magnun.af/src! Aborting..."
    exit
else
    echo "src/ is clean, building..."
fi

cp main_sp4.cpp ../../src
../build.sh
rm ../../src/main_sp4.cpp

if [ -e ../../bin/magnum.af-opencl ];then
    ../../bin/magnum.af-opencl $1
elif [ -e ../../bin/magnum.af-cuda ];then
    ../../bin/magnum.af-cuda $1
else
    ../../bin/magnum.af-cpu $1
fi

if [ -n "$1" ]; then
    ./plot_sp4.sh $1/m.dat
else 
    ./plot_sp4.sh ../../Data/Testing/m.dat
fi
