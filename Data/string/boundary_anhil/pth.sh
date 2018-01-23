#!/bin/bash
mydir=$PWD

#Building
cd ~/git/pth-mag/build
#cmake ..
make -j12
cd $mydir

#mkdir run
rm run/*
cd run
# Running pth-mag-opencl
cp ~/git/pth-mag/src/main* .
cp ~/git/pth-mag/bin/pth-mag-opencl .
screen -d -m bash -c "time ./pth-mag-opencl $PWD > cout.dat"
