#!/bin/bash
gnuplot  -e '
        set terminal pdf;
        set output "x_response.pdf";
        set title "x-Response";
	set xlabel "mu_0*H_x [T]";
	set ylabel "Average Magnetizaion";
	plot "m.dat" u ($5*3.14159*4*1e-7):2 w l t "<m_x>";

        set output "full_elliptical_fields.pdf";

        set title "x-Response";
	set xlabel "mu_0*H_x [T]";
	set ylabel "Average Magnetizaion";
	plot 
        "m.dat" u ($5*3.14159*4*1e-7):2 w l t "<m_x>", 
        "m.dat" u ($5*3.14159*4*1e-7):3 w l t "<m_y>",
        "m.dat" u ($5*3.14159*4*1e-7):4 w l t "<m_z>";

        set title "y-Response";
	set xlabel "mu_0*H_y [T]";
	set ylabel "Average Magnetizaion";
	plot 
        "m.dat" u ($6*3.14159*4*1e-7):2 w l t "<m_x>", 
        "m.dat" u ($6*3.14159*4*1e-7):3 w l t "<m_y>",
        "m.dat" u ($6*3.14159*4*1e-7):4 w l t "<m_z>";

        set title "z-Response";
	set xlabel "mu_0*H_z [T]";
	set ylabel "Average Magnetizaion";
	plot 
        "m.dat" u ($7*3.14159*4*1e-7):2 w l t "<m_x>", 
        "m.dat" u ($7*3.14159*4*1e-7):3 w l t "<m_y>",
        "m.dat" u ($7*3.14159*4*1e-7):4 w l t "<m_z>";

        set notitle;
	set xlabel "step";
	set ylabel "Average Magnetizaion";
	plot 
        "m.dat" u 1:2 w l t "<m_x>", 
        "m.dat" u 1:3 w l t "<m_y>",
        "m.dat" u 1:4 w l t "<m_z>";

	set xlabel "step";
	set ylabel "mu_0*H_ext [T]";
	plot 
        "m.dat" u 1:($5*3.14159*4*1e-7) w l t "H_x",
        "m.dat" u 1:($6*3.14159*4*1e-7) w l t "H_y",
        "m.dat" u 1:($7*3.14159*4*1e-7) w l dt 4 t "H_z";
' --persist
