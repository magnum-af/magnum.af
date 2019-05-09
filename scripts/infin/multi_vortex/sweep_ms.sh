#!/bin/bash
#$1 path to write to
cd "$(dirname "${BASH_SOURCE[0]}")"

nx_disk=100
a=1700
iterategpu=0
for ms in {1500..2000..50}; do # a is the spacing in nm between disks
    echo "ms = $ms"
    selectgpu=$(($iterategpu % 4))
    echo "iterate gpu = $iterategpu"
    echo "select gpu  = $selectgpu"
    ../../magnum.af -S -g "$selectgpu" multi_vortex_demag.py "$1"/"$ms" "$nx_disk" "$a" "$ms"
    ((iterategpu++))
done

cp sweep_ms_postpr.sh "$1"/
