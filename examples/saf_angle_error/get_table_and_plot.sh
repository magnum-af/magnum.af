#!/bin/bash
# write numerically sorted folders' table.dat values to table_combined.dat
ls -v RKKY*/table.dat | xargs cat | tee table_combined.dat
gnuplot plot_table_heatmap.gpi
gnuplot plot_ref_layer_over_files.gpi
