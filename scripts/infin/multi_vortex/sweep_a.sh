#!/bin/bash
#$1 path to write to (ending with "/")
nx_disk=100
iterategpu=0
for a in {500..2000..100}; do # a is the spacing in nm between disks
    echo $a
    selectgpu=$(($iterategpu % 4))
    echo "iterate gpu = $iterategpu"
    echo "select gpu  = $selectgpu"
    ../../magnum.af -S -g "$selectgpu" multi_vortex_demag.py "$1"/"$a" "$nx_disk" "$a"
    ((iterategpu++))
done

cp sweep_a_postpr.sh "$1"/
sed -i "/gnuplot/ a \ \ \ \ set title \"Hx_{demag} as function of spacing a.\";" "$1"/sweep_a_postpr.sh
