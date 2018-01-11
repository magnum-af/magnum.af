#!/bin/bash
gnuplot -e '
 set logscale x;
 set format x "%.1e"; 
 p  "~/git/pth-mag/Data/rigid/rigid.dat" u 1:3  w lp t "T=10K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:4  w lp t "abs T=10K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:5  w l  t "ref T=10K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:7  w lp t "T=50K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:8  w lp t "abs T=50K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:9  w l  t "ref T=50K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:11 w lp t "T=200K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:12 w lp t "abs T=200K",
    "~/git/pth-mag/Data/rigid/rigid.dat" u 1:13 w l  t "ref T=200K"' -persist
