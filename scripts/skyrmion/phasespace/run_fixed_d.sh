#!/bin/bash
# .) compile by executing run_single.sh
# .) cp binary in data/ folder and cd data/
# .) enter screen
# .) run this script

echo "Enter GPU"
read GPU
echo "$GPU"
echo "Enter D_DMI"
read D
echo "$D"
K=10400000
Kstep=800000
mkdir "d$D"
for value in {0..9..1} 
do
  mkdir "d$D/k$K"
  ./magnum.af-cuda "d$D/k$K" "$GPU" "$D" "$K" > d$D/k$K/cout.dat
  K=$(awk "BEGIN {print $K-$Kstep; exit}")
done
