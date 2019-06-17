#!/bin/bash -e
# $1 ... path to magnum.af root directory
# $2 ... verbose. If "$2" != "true", stdout is redirected to /dev/null

currdir=$PWD
if [ -d $1/bin ]
then
    if [ -z "(ls -A $1/bin)" ]; then
        [ "$2" == "true" ] && echo "removing files in $1/bin/*"
        rm $1/bin/*
    fi
fi
if [ ! -d $1/build ]
then
    mkdir $1/build
fi
cd $1/build

[ "$2" == "true" ] && cmake .. || cmake .. > /dev/null

[ "$2" == "true" ] && make -j || make -j > /dev/null
cd $currdir
