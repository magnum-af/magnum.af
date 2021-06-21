#!/bin/bash -ex
# $1 scriptname
# $2 optional outdir
script="$1"
outdir=${2:-outdir_screen_$1/}
mkdir -p $outdir
screen -d -m bash -c "$script -o $outdir |& tee $outdir/cout.log"
