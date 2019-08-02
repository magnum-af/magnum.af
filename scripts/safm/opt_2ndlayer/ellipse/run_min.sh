#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

GPU="1"

# geometry
nx="128"
ny="32"
x="100e-9"
y="25e-9"

# hys info
hzee_max="2.0" # maximum H_external for hysteresis loop
steps="200"

# write dirs
rootdir="$HOME"/data_magnum.af/safm/opt_2ndlayer/ellipse_"$nx"nx_"$ny"ny_"$x"x_"$y"y/
writedir_hf="$rootdir"h_free_layer/
writedir_hys="$rootdir"rk_"$hzee_max"T_"$integr_time"s/

# calculating h_free_layer.vti for mesh if not existant
if [ ! -f "$writedir_hf/h_free_layer.vti" ]; then
    ../../../magnum.af safm_ellipse.cpp "$writedir_hf" "$nx" "$ny" "$x" "$y"
fi


# calculating hysteresis

../../../magnum.af -g "$GPU" -p plot.gpi safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer.vti
