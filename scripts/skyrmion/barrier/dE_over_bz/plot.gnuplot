set terminal pdfcairo enhanced
set xlabel "b_z [J_{ex}]"
set ylabel "dE [J_{ex}]"
set grid
set output "dE_over_bz.pdf"
plot "bz_e_barriers.dat" u 1:($2/(2*1.0e-9*1.6e-11)) w lp notitle
