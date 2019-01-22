#!/bin/bash
# this scripts tests wether there is a main*.cpp file in /src
# if so, it is moved to magnum.af/
# $1 ... verbose

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..


if [ -e $magafdir/src/main*.cpp ]
then
    [ "$1" == "true" ] && echo "Temoraryly moving current main in /src to /temp_main"

    mkdir --parents $magafdir/temp_main
    mv $magafdir/src/main*.cpp $magafdir/temp_main
else
    [ "$1" == "true" ] && echo "src/ is clean, building..."
fi
