#!/bin/bash
# $1 ... path to magnum.af root directory
# $2 ... verbose. If "$2" != "true", stdout is redirected to /dev/null

currdir=$PWD
if [ -d $1/bin ]
then 
    [ "$2" == "true" ] && echo "removing files in $1/bin/*"
    rm $1/bin/*
fi
if [ ! -d $1/build ]
then 
    mkdir $1/build
fi
cd $1/build

if [ "$HOSTNAME" = SERVERGPU1 ]; then
    if [ "$2" == "true" ]; then 
        /home/paul/Programs/cmake-3.13.0-rc2-Linux-x86_64/bin/cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build -DArrayFire_DIR=/usr/local/arrayfire/share/ArrayFire/cmake ..
    else
        /home/paul/Programs/cmake-3.13.0-rc2-Linux-x86_64/bin/cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build -DArrayFire_DIR=/usr/local/arrayfire/share/ArrayFire/cmake .. > /dev/null
    fi
else
    [ "$2" == "true" ] && cmake .. || cmake .. > /dev/null
fi

[ "$2" == "true" ] && make -j || make -j > /dev/null
cd $currdir
