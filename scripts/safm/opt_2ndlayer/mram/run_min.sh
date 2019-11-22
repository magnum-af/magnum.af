#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

GPU="0"

# geometry
x="60e-9"
y="60e-9"
nx="64"
ny="64"

# hys info
hzee_max="1.0" # maximum H_external for hysteresis loop
steps="200"

# write dirs
rootdir="$HOME"/data_magnum.af/safm/opt_2ndlayer/mram/"$x"x_"$y"y_"$nx"nx_"$ny"ny/
writedir_hf="$rootdir"h_free_layer/
writedir_hys="$rootdir"min_"$hzee_max"T_"$steps"steps/

# calculating h_free_layer.vti for mesh if not existant
if [ ! -f "$writedir_hf/h_free_layer.vti" ]; then
    ../../../magnum.af safm_ellipse.cpp "$writedir_hf" "$nx" "$ny" "$x" "$y"
fi


# calculating hysteresis
#../../../magnum.af -g "$GPU" -p plot.gpi safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer.vti

for it in {0..16..2}
do
    writedir_hys="$rootdir"min_"$hzee_max"T_"$steps"steps_it_5degree/it"$it"/
    ../../../magnum.af -g "$GPU" -p plot.gpi safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer_it_"$it".vti
done
