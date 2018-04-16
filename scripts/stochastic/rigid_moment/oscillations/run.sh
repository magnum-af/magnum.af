#!/bin/bash
path=$(pwd)
cd ~/git/pth-mag/build
cmake ..
make -j
cd $path
cp ~/git/pth-mag/src/main* .
cp ~/git/pth-mag/bin/pth-mag-cpu .
time ./pth-mag-cpu $PWD > cout.dat
