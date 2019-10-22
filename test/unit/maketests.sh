#!/bin/bash -e
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # call this scripts directory

if [ -d ./bin ]
then
    if [ -z "(ls -A ./bin)" ]; then
        echo "removing files in ./bin/*"
        rm ./bin/*
    fi
fi
if [ -d ./build ]
then
    if [ -z "(ls -A ./build)" ]; then
        echo "cleaning up existing build/ by removing all files in ./build/*"
        rm ./build/*
    fi
else
    mkdir ./build
fi
cd ./build

if [ "$HOSTNAME" = SERVERGPU1 ]; then
    /home/paul/Programs/cmake-3.13.0-rc2-Linux-x86_64/bin/cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build -DArrayFire_DIR=/usr/local/arrayfire/share/ArrayFire/cmake ..
else
    cmake ..
fi

make -j
cd $currdir
