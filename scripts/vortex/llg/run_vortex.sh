#!/bin/bash
# magnum.af runfile
# Usage:
# ./magnum.af $1=main.cpp $2=absolute/new/output/path $3=optional_GPU-number $4<A>

# Exit on error
set -e

# Setting (input) variables
magafdir=$( dirname "${BASH_SOURCE[0]}" )/../../.. # path to magnum.af/
buildfile=$1
outputpath=$2
if [[ "$3" == "" ]]; then
    GPU=0
else
    GPU=$3
fi
# if given, copying plotfile in outputdir
if [[ "$buildfile" == "" ]]; then
    echo "Error: missing buildfile as first argument required, aborting..."
    exit 1
fi
# checking if other main exists in /src
$magafdir/scripts/bash/check_main.sh

# building
cp $buildfile $magafdir/src
$magafdir/scripts/bash/build.sh $magafdir
rm $magafdir/src/$(basename $buildfile)

# moving possible old main back
$magafdir/scripts/bash/check_main_remove.sh


for valueAB in {0..100..10}
do
    # copying files
    mkdir -p $outputpath/$valueAB
    if [[ "$outputpath/$valueAB" != "" ]]; then
        cp $magafdir/bin/magnum.af-* $outputpath/$valueAB
        cp $buildfile $outputpath/$valueAB
        cp plot_hysteresis.sh $outputpath/$valueAB
    fi
    # creating write dir
    echo $valueAB
    $magafdir/bin/magnum.af-cuda $outputpath/$valueAB $GPU $4 $valueAB > $outputpath/$valueAB/cout.txt 2>&1
    #$magafdir/bin/magnum.af-cpu $outputpath/$valueAB $GPU $4 $valueAB > $outputpath/$valueAB/cout.txt 2>&1
done
