#!/bin/bash
# this scripts tests wether there is a main*.cpp file in /scripts/single_script_build
# if so, it is moved to magnum.af/
# $1 ... verbose

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..


if [ -e $magafdir/scripts/single_script_build/main*.cpp ]
then
    [ "$1" == "true" ] && echo "Temoraryly moving current main in /scripts/single_script_build to /temp_main"

    mkdir --parents $magafdir/temp_main
    mv $magafdir/scripts/single_script_build/main*.cpp $magafdir/temp_main
else
    [ "$1" == "true" ] && echo "scripts/single_script_build/ is clean, building..."
fi
exit 0
