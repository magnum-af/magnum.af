#!/bin/bash
simtime=100 # in [ns]
nstep=1000 # mean every nth step
../magnum.af -fs -p plot.gnu offset_sensor.py output "$simtime" "$nstep"
