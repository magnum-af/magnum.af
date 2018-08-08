#!/bin/bash
# Usage .sh absolute/path/to/write/output/ <optional-GPU number> /path/to/mrelax.vti
set -e

# calling this scripts's directory
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# relative path to magnum.af/
magafdir=../..
buildfile=main_elliptical_minimizer.cpp
plotfile=../vortex/plot_hysteresis.sh

# checking if other main exists in /src
$magafdir/scripts/bash/check_main.sh

# building 
cp $buildfile $magafdir/src
$magafdir/scripts/bash/build.sh $magafdir
rm $magafdir/src/$buildfile

# moving possible old main back
$magafdir/scripts/bash/check_main_remove.sh

# creating write dir
$magafdir/scripts/bash/check_write_dir.sh $1
if [ "$?" == 1 ];then
    exit 1
fi

# copying files
cp $magafdir/bin/magnum.af-* $1
cp $buildfile $1
cp $plotfile $1

# running
if [ -e $magafdir/bin/magnum.af-cuda ];then
    screen -d -S GPU$2 -m bash -c "export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64:$LD_LIBRARY_PATH && export PATH=/usr/local/cuda-8.0/bin:$PATH && $magafdir/bin/magnum.af-cuda $1 $2 $3> $1/cout.txt 2>&1"
elif [ -e $magafdir/bin/magnum.af-opencl ];then
    screen -d -S GPU$2 -m bash -c "$magafdir/bin/magnum.af-opencl $1 $2 $3 > $1/cout.txt 2>&1"
    echo "NOTE: Opencl Backend fails due to arrayfire error 202 in fft, if chosen mesh dimensions have primefactor > 11"
else
    $magafdir/bin/magnum.af-cpu $1 $2 $3
fi

echo "To follow cout.dat run:"
echo "tail -f $1/cout.txt"
