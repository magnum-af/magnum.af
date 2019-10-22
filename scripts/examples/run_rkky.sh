#!/bin/bash -e

magafdir=../../
outputdir=${1:-"output_rkky/"}

# check if binaries exist and build otherwise
for f in $magafdir/bin/rkky_*; do
    [ -e "$f" ] && echo "binaries do exist" || ( echo "binaries do not exist" && (cd $magafdir/build && cmake .. && make ) )
    break
done

# mkdir -p
[ ! -d "$outputdir" ] && mkdir -p $outputdir

# set binary, preferring cuda over opencl over cpu
binary="$magafdir"/bin/rkky_cpu
[ -e "$magafdir"/bin/rkky_opencl ] && binary="$magafdir"/bin/rkky_opencl
[ -e "$magafdir"/bin/rkky_cuda ] && binary="$magafdir"/bin/rkky_cuda

# run
$binary $outputdir

# plot
gnuplot -e "
    p '$outputdir/E.dat';
    set title sprintf('E over angle, dE=%.3e', (GPVAL_DATA_Y_MAX-GPVAL_DATA_Y_MIN));
    set terminal pdfcairo;
    set output '$outputdir/E.pdf';
    p '$outputdir/E.dat' u 1:2 w lp t 'E', GPVAL_DATA_Y_MAX t sprintf('E_{max}=%.3e', GPVAL_DATA_Y_MAX), GPVAL_DATA_Y_MIN t sprintf('E_{min}=%.3e', GPVAL_DATA_Y_MIN);
    set terminal pop;
    replot;
" --persist
