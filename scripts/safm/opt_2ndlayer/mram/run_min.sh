#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

GPU="3"

# geometry
x="60e-9"
y="60e-9"
nx="64"
ny="64"

# hys info
hzee_max="1.0" # maximum H_external for hysteresis loop
steps="2000"

# write dirs
rootdir="$HOME"/data_magnum.af/safm/opt_2ndlayer/mram/eval_demag_only_in_material_init5.001_hfl_only_in_mat/"$x"x_"$y"y_"$nx"nx_"$ny"ny/
writedir_hf="$rootdir"h_free_layer/
writedir_hys="$rootdir"min_"$hzee_max"T_"$steps"steps/

# calculating h_free_layer.vti for mesh if not existant
if [ ! -f "$writedir_hf/h_free_layer.vti" ]; then
    ../../../magnum.af safm_ellipse.cpp "$writedir_hf" "$nx" "$ny" "$x" "$y"
fi


# calculating hysteresis
#../../../magnum.af -g "$GPU" -p plot.gpi safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer.vti

writedir="$rootdir"min_"$hzee_max"T_"$steps"steps_it_5degree_zee/
mkdir "$writedir"
cp "${BASH_SOURCE[0]}" "$writedir"
cp plot_multi_hys.gpi "$writedir"
for it in {0..14..2} # specific range for one example
do
    writedir_hys="$writedir"it"$it"/
    #../../../magnum.af -g "$GPU" -p plot.gpi safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer_it_"$it".vti
    ../../../magnum.af -g "$GPU" safm_ellipse_hys_minimizer.cpp "$writedir_hys" "$hzee_max" "$steps" "$writedir_hf"/h_free_layer_it_"$it".vti
done
cd "$writedir"
gnuplot plot_multi_hys.gpi
