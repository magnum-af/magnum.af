#!/bin/bash
#$1 path to write to
cd "$(dirname "${BASH_SOURCE[0]}")"

nx_disk=100
a=1700
iterategpu=0
ms=1750
for degree in {0..90..10}; do # a is the spacing in nm between disks
    echo "degree = $degree"
    selectgpu=$(($iterategpu % 4))
    echo "iterate gpu = $iterategpu"
    echo "select gpu  = $selectgpu"
    ../../magnum.af -S -g "$selectgpu" multi_vortex_demag.py "$1"/"$degree" "$nx_disk" "$a" "$ms" "$degree"
    ((iterategpu++))
done

cp sweep_angle_postpr.sh "$1"/
sed -i "/gnuplot/ a \ \ \ \ set title \"Angular dependency of in-plain magnetized disks. \\\n Ms*mu_0=$ms [mT], a=$a [nm]\";" "$1"/sweep_angle_postpr.sh
