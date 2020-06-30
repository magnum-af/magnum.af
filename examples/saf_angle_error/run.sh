#!/bin/bash

(cd ../../build/ && make -j)
Hmax_mT="0.100"
Jaf="0.36e-3"
#Jaf="0.49e-3"
#Jaf="0.72e-3"
#Jaf="1.08e-3"

Ms1_start=500000
Ms1_step=100000
Ms1_stop=2000000
RKKY_start=-0
RKKY_step=-0.1
RKKY_stop=-1.6
parent_outpath=$HOME/data_magnum.af/inf_saf/"Jaf$Jaf"/

for Ms1 in $(LC_ALL=C seq "$Ms1_start" "$Ms1_step" "$Ms1_stop" ); do
    for RKKY in $(LC_ALL=C seq "$RKKY_start" "$RKKY_step" "$RKKY_stop"); do
        echo "$Ms1, $RKKY, $Jaf"
        outpath="$parent_outpath"RKKYmT"$RKKY"Ms1"$Ms1"/
        echo "$outpath"
        mkdir -p "$outpath"
        GPU="0"
        ../../bin/saf_angle_error "$outpath" "$GPU" "$Hmax_mT" "$RKKY" "$Ms1" "$Jaf" | tee "$outpath"log.txt
    done
done

# postprocessing:
cp plot_ref_layer_over_files.gpi "$parent_outpath"
cp plot_table_heatmap.gpi "$parent_outpath"
cp get_table_and_plot.sh "$parent_outpath"
cd "$parent_outpath"
./get_table_and_plot.sh
