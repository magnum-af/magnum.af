#!/bin/bash
cp ~/git/magnum.af/src/main*.cpp .
cp ~/git/magnum.af/bin/magnum.af-opencl .
mkdir data
screen -d -m bash -c 'time ./magnum.af-opencl $PWD/data 2 > $PWD/data/cout.dat'
