#!/bin/bash
gnuplot -e '
 set logscale x;
 set format x "%.1e"; 
 set xlabel "dt";
 set yrange [-0.6:1.1];
 set ylabel "<|m_z|>";
 p  
    "rigid_Heun/rigid.dat"         u 1:4  w p lt 1  dashtype 1 t "    T=10K",
    "rigid_Heun/rigid.dat"         u 1:5  w l lt 1  dashtype 3 t "ref T=10K",
    "rigid_Heun/rigid.dat"         u 1:8  w p lt 2  dashtype 1 t "    T=50K",
    "rigid_Heun/rigid.dat"         u 1:9  w l lt 2  dashtype 3 t "ref T=50K",
    "rigid_Heun/rigid.dat"         u 1:12 w p lt 3  dashtype 1 t "    T=200K",
    "rigid_Heun/rigid.dat"         u 1:13 w l lt 3  dashtype 3 t "ref T=200K",
    "rigid_Heun_cpu/rigid.dat"     u 1:4  w p lt 4  dashtype 1 t "    T=10K",
    "rigid_Heun_cpu/rigid.dat"     u 1:5  w l lt 4  dashtype 3 t "ref T=10K",
    "rigid_Heun_cpu/rigid.dat"     u 1:8  w p lt 5  dashtype 1 t "    T=50K",
    "rigid_Heun_cpu/rigid.dat"     u 1:9  w l lt 5  dashtype 3 t "ref T=50K",
    "rigid_Heun_cpu/rigid.dat"     u 1:12 w p lt 6  dashtype 1 t "    T=200K",
    "rigid_Heun_cpu/rigid.dat"     u 1:13 w l lt 6  dashtype 3 t "ref T=200K",
    "rigid_SemiHeun/rigid.dat"     u 1:4  w p lt 7  dashtype 1 t "    T=10K",
    "rigid_SemiHeun/rigid.dat"     u 1:5  w l lt 7  dashtype 3 t "ref T=10K",
    "rigid_SemiHeun/rigid.dat"     u 1:8  w p lt 8  dashtype 1 t "    T=50K",
    "rigid_SemiHeun/rigid.dat"     u 1:9  w l lt 8  dashtype 3 t "ref T=50K",
    "rigid_SemiHeun/rigid.dat"     u 1:12 w p lt 9  dashtype 1 t "    T=200K",
    "rigid_SemiHeun/rigid.dat"     u 1:13 w l lt 9  dashtype 3 t "ref T=200K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:4  w p lt 10 dashtype 1 t "    T=10K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:5  w l lt 10 dashtype 3 t "ref T=10K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:8  w p lt 11 dashtype 1 t "    T=50K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:9  w l lt 11 dashtype 3 t "ref T=50K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:12 w p lt 12 dashtype 1 t "    T=200K",
    "rigid_SemiHeun_cpu/rigid.dat" u 1:13 w l lt 12 dashtype 3 t "ref T=200K"
' -persist

#gnuplot -e '
# set logscale x;
# set format x "%.1e"; 
# set xlabel "dt";
# set yrange [-0.6:1.1];
# set ylabel "<|m_z|>";
# p  
#    "rigid_Heun/rigid.dat"         u 1:4  w p lt 1  dashtype 1 t "    T=10K",
#    "rigid_Heun/rigid.dat"         u 1:5  w l lt 1  dashtype 3 t "ref T=10K",
#    "rigid_Heun/rigid.dat"         u 1:8  w p lt 2  dashtype 1 t "    T=50K",
#    "rigid_Heun/rigid.dat"         u 1:9  w l lt 2  dashtype 3 t "ref T=50K",
#    "rigid_Heun/rigid.dat"         u 1:12 w p lt 3  dashtype 1 t "    T=200K",
#    "rigid_Heun/rigid.dat"         u 1:13 w l lt 3  dashtype 3 t "ref T=200K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:4  w p lt 4  dashtype 1 t "    T=10K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:5  w l lt 4  dashtype 3 t "ref T=10K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:8  w p lt 5  dashtype 1 t "    T=50K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:9  w l lt 5  dashtype 3 t "ref T=50K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:12 w p lt 6  dashtype 1 t "    T=200K",
#    "rigid_Heun_cpu/rigid.dat"     u 1:13 w l lt 6  dashtype 3 t "ref T=200K",
#' -persist

#gnuplot -e '
# set terminal pdfcairo size 8,5.5;
# set output "rigid.pdf";
# set logscale x;
# set format x "%.1e"; 
# set xlabel "dt";
# set ylabel "<|m_z|>";
# set yrange [-0.6:1.1];
# p  
#    "rigid.dat" u 1:4  w p lt 1 lw 1 t "    T=10K",
#    "rigid.dat" u 1:5  w l lt 1 lw 2 t "ref T=10K",
#    "rigid.dat" u 1:8  w p lt 2 lw 1 t "    T=50K",
#    "rigid.dat" u 1:9  w l lt 2 lw 2 t "ref T=50K",
#    "rigid.dat" u 1:12 w p lt 3 lw 1 t "    T=200K",
#    "rigid.dat" u 1:13 w l lt 3 lw 2 t "ref T=200K"'