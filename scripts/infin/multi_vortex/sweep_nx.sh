#!/bin/bash

for i in {50..250..10}; do
    echo $i
    /home/paul/magnum.af/scripts/magnum.af multi_vortex_demag.py sweep_nx/$i $i
done
