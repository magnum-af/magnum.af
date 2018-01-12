#!/bin/bash
gnuplot -e '
 set logscale x;
 set format x "%.1e"; 
 set xlabel "dt";
 set yrange [-0.6:1.1];
 set ylabel "<|m_z|>";
 p  
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:4  w p lt 1 dashtype 1 t "    T=10K",
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:5  w l lt 1 dashtype 3 t "ref T=10K",
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:8  w p lt 2 dashtype 1 t "    T=50K",
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:9  w l lt 2 dashtype 3 t "ref T=50K",
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:12 w p lt 3 dashtype 1 t "    T=200K",
    "~/git/pth-mag/Data/randu_rigid/rigid.dat" u 1:13 w l lt 3 dashtype 3 t "ref T=200K"' -persist

#gnuplot -e '
# set terminal pdfcairo size 8,5.5;
# set output "~/git/pth-mag/Data/rigid/rigid.pdf";
# set logscale x;
# set format x "%.1e"; 
# set xlabel "dt";
# set ylabel "<|m_z|>";
# set yrange [-0.6:1.1];
# p  
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:4  w p lt 1 dashtype 1 lw 1 t "    T=10K",
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:5  w l lt 1 dashtype 3 lw 2 t "ref T=10K",
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:8  w p lt 2 dashtype 1 lw 1 t "    T=50K",
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:9  w l lt 2 dashtype 3 lw 2 t "ref T=50K",
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:12 w p lt 3 dashtype 1 lw 1 t "    T=200K",
#    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:13 w l lt 3 dashtype 3 lw 2 t "ref T=200K"'
#




    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:3  w lp t "T=10K",
    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:7  w lp t "T=50K",
    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:11 w lp t "T=200K",
    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:3  w lp t "T=10K",
    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:7  w lp t "T=50K",
    #"~/git/pth-mag/Data/rigid/rigid.dat" u 1:11 w lp t "T=200K",
