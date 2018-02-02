#!/bin/bash
#echo $PWD
for dir in ./*/
do
    cd ${dir}
    #echo $PWD
    sbatch ./vsc-parallel-pth-mag-cpu.sh
    cd ..
done
