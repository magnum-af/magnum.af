#!/bin/bash

(cd ../../build/ && make -j)
GPU="0"
Hmax_mT="0.100"
#Jaf="0.36e-3"
#Jaf="0.54e-3"
#Jaf="0.72e-3"
#Jaf="0.90e-3"
#Jaf="1.08e-3"

Ms1_start=500000
Ms1_step=100000
Ms1_stop=2000000
RKKY_start=-0
RKKY_step=-0.1
RKKY_stop=-1.6
parent_outpath=$HOME/data_magnum.af/inf_saf/"Jaf$Jaf"/

for Ms1 in $(LC_ALL=C seq "$Ms1_start" "$Ms1_step" "$Ms1_stop" ); do
    for RKKY in $(LC_ALL=C seq "$RKKY_start" "$RKKY_step" "$RKKY_stop"); do
        echo "$Ms1, $RKKY, $Jaf"
        outpath="$parent_outpath"RKKYmT"$RKKY"Ms1"$Ms1"/
        echo "$outpath"
        if [ ! -d "$outpath" ]; then
            mkdir -p "$outpath"
            time ../../bin/saf_angle_error "$outpath" "$GPU" "$Hmax_mT" "$RKKY" "$Ms1" "$Jaf" | tee "$outpath"log.txt
        fi
    done
done

# postprocessing:
cp plot.gpi "$parent_outpath"
cd "$parent_outpath"
# write numerically sorted folders' table.dat values to table_combined.dat
ls -v RKKY*/table.dat | xargs cat | tee table_combined.dat
gnuplot plot.gpi

# getting min of mean(abs(my)) values
min=$(sed '/^#/d' < table_combined.dat | awk '{print $4}' | sort -g | head -1)
cat table_combined.dat | grep "$min" | tee min_mean_abs_value.dat

# plot min(mean(abs(my))) over J_af:
cd ..
cat Jaf*/min_mean_abs_value.dat > minvals.dat

echo "set terminal pdf
set output 'fig_min_over_Jaf.pdf'
set xlabel 'J_{af} [J/m2]'
set ylabel 'min(mean(abs(m_y)))'
p 'minvals.dat' u 5:4 w lp notitle
set terminal jpeg
set output 'fig_min_over_Jaf.jpg'
replot" > plot_min_over_Jaf.gpi

gnuplot plot_min_over_Jaf.gpi
