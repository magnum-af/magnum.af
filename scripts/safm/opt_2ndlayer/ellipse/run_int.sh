#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"

writedir="$HOME"/data_magnum.af/safm/opt_2nd/ellipse
hzee_max="0.08" # maximum H_external for hysteresis loop
#integr_time="200e-9"
integr_time="800e-9"

if [ ! -f "$writedir/h_free_layer.vti" ]; then
    ../../../magnum.af -f safm_ellipse.cpp "$writedir"
fi

../../../magnum.af -fp plot.gpi safm_ellipse_hys_rk.cpp "$writedir" "$hzee_max" "$integr_time"
