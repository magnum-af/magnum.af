#!/bin/bash
#gnuplot -e 'p "m.dat" u 1:2 w l t "mx","m.dat" u 1:3 w l t "my", "m.dat" u 1:4 w l t "mz"' -persist
#gnuplot -e 'p "m.dat" u 1:2 w l t "mx","m.dat" u 1:3 w l t "my", "m.dat" u 1:4 w l t "mz","temp_m.dat" u 1:2 w l t "mx","temp_m.dat" u 1:3 w l t "my", "temp_m.dat" u 1:4 w l t "mz"' -persist
#gnuplot -e 'p "m.dat" u 1:2 w l t "mx","m.dat" u 1:3 w l t "my", "m.dat" u 1:4 w l t "mz","temp_m.dat" u 1:2 w l t "mx","temp_m.dat" u 1:3 w l t "my", "temp_m.dat" u 1:4 w l t "mz","rk4_m.dat" u 1:2 w l t "mx","rk4_m.dat" u 1:3 w l t "my", "rk4_m.dat" u 1:4 w l t "mz"' -persist
#gnuplot -e 'p "m.dat" u 1:2 w l t "mx","m.dat" u 1:3 w l t "my", "m.dat" u 1:4 w l t "mz","rk4_m.dat" u 1:2 w l t "mx","rk4_m.dat" u 1:3 w l t "my", "rk4_m.dat" u 1:4 w l t "mz"' -persist
gnuplot -e 'p "~/git/magnum.af/Data/Testing/m.dat" u 1:2 w l t "mx","~/git/magnum.af/Data/Testing/m.dat" u 1:3 w l t "my", "~/git/magnum.af/Data/Testing/m.dat" u 1:4 w l t "mz","~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:2 w l t "mx","~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:3 w l t "my", "~/git/magnum.af/Data/Testing/rk4_m.dat" u 1:4 w l t "mz"' -persist
