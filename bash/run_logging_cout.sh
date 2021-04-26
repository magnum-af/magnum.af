#!/bin/bash -ex
# script logging output from cout to file cout.txt
# $1 executable to run
# $2 output dir to write cout.txt to

script=$1
outdir=${2:-"output_$1"}
mkdir -p "$outdir"
./"$script" -o "$outdir" | tee "$outdir"/cout.txt
