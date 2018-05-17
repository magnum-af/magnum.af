#!/bin/bash
path=$(pwd)
cd ~/git/magnum.af/build
cmake ..
make -j
cd $path
cp ~/git/magnum.af/src/main* .
cp ~/git/magnum.af/bin/magnum.af-cpu .
time ./magnum.af-cpu $PWD > cout.dat
