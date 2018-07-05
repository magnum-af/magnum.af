#!/bin/bash
gnuplot -e 'p "'$1'" u 1:2 w l t "mx","'$1'" u 1:3 w l t "my", "'$1'" u 1:4 w l t "mz"' -persist 

#',"~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:2 w l t "mx","~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:3 w l t "my", "~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:4 w l t "mz"' -persist
