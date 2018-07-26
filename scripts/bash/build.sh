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
make -j
cd $currdir
