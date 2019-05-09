#!/bin/bash
#$1 path to write to
cd "$(dirname "${BASH_SOURCE[0]}")"

nx_disk=100
a_y=1700
ms=1750
degree=0

iterategpu=0
for ax_ay_ratio in {100..300..20}; do # a_x/a_y ratio in %
    echo "ax_ay_ratio = $ax_ay_ratio"
    selectgpu=$(($iterategpu % 4))
    echo "select gpu  = $selectgpu"
    ../../magnum.af -S -g "$selectgpu" multi_vortex_demag.py "$1"/"$ax_ay_ratio" "$nx_disk" "$a_y" "$ms" "$degree" "$ax_ay_ratio"
    ((iterategpu++))
done

cp sweep_ax_postpr.sh "$1"/
sed -i "/gnuplot/ a \ \ \ \ set title \"x-Spacing Dependency. \\\n Ms*mu_0=$ms [mT], a_y=$a_y [nm]\";" "$1"/sweep_ax_postpr.sh
