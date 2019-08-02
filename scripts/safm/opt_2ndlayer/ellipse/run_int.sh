#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

hzee_max="0.25" # maximum H_external for hysteresis loop
integr_time="300e-9"
writedir="$HOME"/data_magnum.af/safm/opt_2ndlayer/ellipse_rk_"$hzee_max"T_"$integr_time"s

if [ ! -f "$writedir/h_free_layer.vti" ]; then
    ../../../magnum.af -f safm_ellipse.cpp "$writedir"
fi

../../../magnum.af -fp plot.gpi safm_ellipse_hys_rk.cpp "$writedir" "$hzee_max" "$integr_time"
