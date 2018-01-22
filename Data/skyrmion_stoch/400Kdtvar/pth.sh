#!/bin/bash
mydir=$PWD
#dt="1e-13"
echo "type dt"
read dt
echo $dt
#GPU=0
echo "type GPU"
read GPU
echo $GPU
T=400
#echo "type T"
#read T
#echo $T

#Building
cd ~/git/pth-mag/build
make -j12
cd $mydir

# Running pth-mag-opencl
mkdir $dt
cd $dt
mkdir skyrm
cp ~/git/pth-mag/src/main_skyrmion_stoch.cpp .
cp ~/git/pth-mag/bin/pth-mag-opencl .
screen -d -m bash -c "time ./pth-mag-opencl $PWD $GPU $T $dt > cout.dat"
