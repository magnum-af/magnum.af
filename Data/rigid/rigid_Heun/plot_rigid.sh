#!/bin/bash
gnuplot -e '
 set logscale x;
 set format x "%.1e"; 
 set xlabel "dt";
 set yrange [-0.6:1.1];
 set ylabel "<|m_z|>";
 p  
    "rigid.dat" u 1:4  w p lt 1 t "    T=10K",
    "rigid.dat" u 1:5  w l lt 1 t "ref T=10K",
    "rigid.dat" u 1:8  w p lt 2 t "    T=50K",
    "rigid.dat" u 1:9  w l lt 2 t "ref T=50K",
    "rigid.dat" u 1:12 w p lt 3 t "    T=200K",
    "rigid.dat" u 1:13 w l lt 3 t "ref T=200K"' -persist

gnuplot -e '
 set terminal pdfcairo size 8,5.5;
 set output "rigid.pdf";
 set logscale x;
 set format x "%.1e"; 
 set xlabel "dt";
 set ylabel "<|m_z|>";
 set yrange [-0.6:1.1];
 p  
    "rigid.dat" u 1:4  w p lt 1 lw 1 t "    T=10K",
    "rigid.dat" u 1:5  w l lt 1 lw 2 t "ref T=10K",
    "rigid.dat" u 1:8  w p lt 2 lw 1 t "    T=50K",
    "rigid.dat" u 1:9  w l lt 2 lw 2 t "ref T=50K",
    "rigid.dat" u 1:12 w p lt 3 lw 1 t "    T=200K",
    "rigid.dat" u 1:13 w l lt 3 lw 2 t "ref T=200K"'
