#!/bin/bash
simtime=100.0 # in [ns]
nstep=1000 # mean every nth step

../magnum.af -sg 3 -p plot.gnu offset_sensor.py "$1" "$simtime" "$nstep"

## or use docker
#../magnum.af.docker -dtg 3 offset_sensor.py "$1" "$simtime" "$nstep"
