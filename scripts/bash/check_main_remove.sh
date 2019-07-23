#!/bin/bash
# this scripts tests whether there is a main*.cpp file in /temp_main
# if so, it is moved to magnum.af/
# $1 ... verbose

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..

for file in $magafdir/scripts/single_script_build/main*.cpp; do
    [ -e "$file" ] && echo "Warning: in check_main_remove.sh: found main*.cpp in magnum.af/scripts/single_script_build but it should be empty. Consider cleaning up mainfiles in magnum.af/scripts/single_script_build."
    break
done

if [ -e $magafdir/temp_main ];then
    for file in $magafdir/temp_main/main*.cpp; do
        [ "$1" == "true" ] && echo "Moving back file $file from temp_main/ to /scripts/single_script_build"
        mv $magafdir/temp_main/main*.cpp $magafdir/scripts/single_script_build
    done
    rmdir $magafdir/temp_main
fi
