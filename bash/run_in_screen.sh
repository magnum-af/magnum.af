#!/bin/bash -ex
# $1 scriptname
# $2 optional outdir
script="$1"
outdir=${2:-outdir_screen_$1/}
mkdir -p $outdir
runcommand="python3 $script -o $outdir ${@:3} |& tee $outdir/cout.log"
screen -d -m bash -c "$runcommand"
echo "tail -f $outdir/cout.log"
