#!/bin/bash
# magnum.af runfile
# Usage:
# ./magnum.af $1=main.cpp $2=absolute/new/output/path $3=optional_GPU-number $4<A>
#e.g
#./run_AB_sweep.sh main_infineon_CoFeB_oop_fl_elliptical_field.cpp /mnt/afa30f69-bd5b-4236-9ff9-1eb489c3fbfd/data/infineon/CoFeB/elliptical/t200ns/A50mT/ 0 50

# Exit on error
set -e

# Setting (input) variables
magafdir=$( dirname "${BASH_SOURCE[0]}" )/../.. # path to magnum.af/
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

start_val=20
stop_val=100
step=20
valueAB=$start_val
while [  $valueAB -le "$stop_val" ]; do
    echo "Creating dirs: $valueAB"
    mkdir -p $outputpath/$valueAB
    cp $magafdir/bin/magnum.af-* $outputpath/$valueAB
    cp $buildfile $outputpath/$valueAB
    let valueAB=valueAB+$step
done

valueAB=$start_val
while [  $valueAB -le "$stop_val" ]; do
    echo "running $valueAB"
    ls $outputpath/$valueAB
    $outputpath/$valueAB/magnum.af-opencl $outputpath/$valueAB $GPU $4 $valueAB > $outputpath/$valueAB/cout.txt 2>&1
    let valueAB=valueAB+$step
done
