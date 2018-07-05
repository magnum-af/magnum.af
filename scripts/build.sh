#!/bin/bash

currdir=$PWD
echo $PWD
if [ ! -f ~/git/magnum.af/build ]
then 
    mkdir ~/git/magnum.af/build
fi
cd ~/git/magnum.af/build
cmake ..
make -j
cd $currdir
