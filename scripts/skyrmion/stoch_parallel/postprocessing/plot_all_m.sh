#!/bin/bash

gnuplot -e "
set terminal pdfcairo;
list = system('ls $1/m*.dat');

set output '$1/all_m_in_one_plot.pdf';
plot for [file in list] file u 0:1 w l t file;

set output '$1/all_m_in_separate_plot.pdf';
do for [file in list]{plot file u 0:1 w l};
"
