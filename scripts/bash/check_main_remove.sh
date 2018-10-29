#!/bin/bash
# this scripts tests whether there is a main*.cpp file in /temp_main
# if so, it is moved to magnum.af/

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..

for file in $magafdir/src/main*.cpp; do
    [ -e "$file" ] && echo "Warning: in check_main_remove.sh: found main*.cpp in magnum.af/src but it should be empty. Consider cleaning up mainfiles in magnum.af/src."
    break
done

if [ -e $magafdir/temp_main ];then
    for file in $magafdir/temp_main/main*.cpp; do
        echo "Info: Moving back file $file from temp_main/ to /src"
        mv $magafdir/temp_main/main*.cpp $magafdir/src
    done
    rmdir $magafdir/temp_main 
fi
