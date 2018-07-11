#!/bin/bash

currdir=$PWD
echo $PWD
if [ ! -f $1/build ]
then 
    mkdir $1/build
fi
cd $1/build
cmake ..
make -j
cd $currdir
