#!/bin/bash
gnuplot  -e '
        set terminal pdf;
	set xlabel "mu_0*H_ext [T]";
	set ylabel "Average Magnetizaion";
        set output "hys.pdf";
	p 
        "m.dat" u ($5*3.14159*4*1e-7):2 w l t "<m_x>", 
        "m.dat" u ($5*3.14159*4*1e-7):3 w l t "<m_y>",
        "m.dat" u ($5*3.14159*4*1e-7):4 w l t "<m_z>"'
gnuplot  -e '
	set xlabel "mu_0*H_ext [T]";
	set ylabel "Average Magnetizaion";
	p 
        "m.dat" u ($5*3.14159*4*1e-7):2 w l t "<m_x>", 
        "m.dat" u ($5*3.14159*4*1e-7):3 w l t "<m_y>",
        "m.dat" u ($5*3.14159*4*1e-7):4 w l t "<m_z>";' --persist
