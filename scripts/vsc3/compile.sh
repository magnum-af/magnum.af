#!/bin/bash
currdir=$PWD
source /home/lv70895/abert3v/magnumfe.custom
module load cuda/8.0.61
cd ~/git/pth-mag/build
cmake -DArrayFire_DIR=/home/lv70895/heistracher/programs/arrayfire/share/ArrayFire/cmake ..
make -j
cd $currdir

