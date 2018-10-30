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

if [ "$HOSTNAME" = SERVERGPU1 ]; then
    /home/paul/Programs/cmake-3.13.0-rc2-Linux-x86_64/bin/cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build -DArrayFire_DIR=/usr/local/arrayfire/share/ArrayFire/cmake ..
else
    cmake ..
fi

make -j
cd $currdir
