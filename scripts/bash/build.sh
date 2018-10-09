#!/bin/bash

currdir=$PWD
if [ -d $1/bin ]
then 
    echo "removing files in $1/bin/*"
    rm $1/bin/*
fi
if [ ! -f $1/build ]
then 
    mkdir $1/build
fi
cd $1/build
cmake ..
# NOTE: GTO use:
# cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build ..

make -j
cd $currdir
