set terminal pdfcairo enhanced
set output "mx_over_t.pdf"
set xlabel "t [s]"
set ylabel "<mx>"
p "m.dat" u 1:2 w lp notitle
