#!/bin/bash
#echo $PWD
for dir in ./*/
do
    ./postprocessing/run.sh ${dir}/data/anihilationtime.dat ${dir}
done
