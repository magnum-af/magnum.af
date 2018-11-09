#!/bin/bash
# magnum.af runfile
# Usage:
# ./magnum.af $1=main.cpp $2=absolute/new/output/path $3=optional_GPU-number $4<A> 
# Example
# ./magnum.af main.cpp $PWD/run1 0 plotfile 

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


# copying plotfile if existing
#if [[ "$4" != "" ]]; then
#    cp $4 $outputpath
#fi

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
    #if [[ "$outputpath/$valueAB" != "" ]]; then
    #    cp $magafdir/bin/magnum.af-* $outputpath/$valueAB
    #    cp $buildfile $outputpath/$valueAB
    #fi
    # creating write dir
    echo $valueAB
    $magafdir/bin/magnum.af-cuda $outputpath/$valueAB $GPU $4 $valueAB > $outputpath/$valueAB/cout.txt 2>&1
done
## running
#if [ -e $magafdir/bin/magnum.af-cuda ];then
#    echo "starting magnum.af-cuda in screen."
#    echo "To follow cout.dat run:"
#    echo "tail -f $outputpath/cout.txt"
#    screen -d -S GPU$GPU -m bash -c "export LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64:$LD_LIBRARY_PATH && export PATH=/usr/local/cuda-9.0/bin:$PATH && $magafdir/bin/magnum.af-cuda $outputpath $GPU $5 > $outputpath/cout.txt 2>&1"
#    screen -ls
#elif [ -e $magafdir/bin/magnum.af-opencl ];then
#    echo "starting magnum.af-opencl in screen."
#    echo "To follow cout.dat run:"
#    echo "tail -f $outputpath/cout.txt"
#    screen -d -S GPU$GPU -m bash -c "$magafdir/bin/magnum.af-opencl $outputpath $GPU $5 > $outputpath/cout.txt 2>&1"
#    screen -ls
#else
#    echo "starting magnum.af-cpu."
#    $magafdir/bin/magnum.af-cpu $outputpath $GPU $5
#fi

