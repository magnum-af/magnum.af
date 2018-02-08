#!/bin/bash
cp ~/git/pth-mag/src/main*.cpp .
cp ~/git/pth-mag/bin/pth-mag-opencl .
mkdir data
screen -d -m bash -c 'time ./pth-mag-opencl $PWD/data 2 > $PWD/data/cout.dat'
