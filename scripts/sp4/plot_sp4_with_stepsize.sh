#/bin/bash
gnuplot -e '
    set ytics nomirror; 
    set y2tics;
    set xlabel "time [s]";
    set ylabel "m components";
    set y2label "stepsize [s]";
    Dx(y) = ($0==0)?(Dy1=y, 1/0) : (Dy2 = Dy1, Dy1 = y, (Dy1-Dy2)); 
    p "m.dat" u 1:2 w l t "mx", "m.dat" u 1:3 w l t "my", "m.dat" u 1:4 w l t "mz", "m.dat"  using 1:(Dx($1)) axis x1y2 w l t "stepsize";
' -persist
