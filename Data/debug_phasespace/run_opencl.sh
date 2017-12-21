#!/bin/bash
echo "Enter GPU"
read GPU
echo "$GPU"
echo "Enter D_DMI"
read D
echo "$D"
#D=0.01728
#K=6400000
K=9600000
Kstep=1600000
mkdir "d$D"
for value in {0..4..1} 
do
  mkdir "d$D/k$K"
  time ./pth-mag-opencl "d$D/k$K" "$GPU" "$D" "$K" > d$D/k$K/cout.dat
  #time ./pth-mag-cuda "d$D/k$K" "$GPU" "$D" "$K" > d$D/k$K/cout.dat
  K=$(awk "BEGIN {print $K-$Kstep; exit}")
done
