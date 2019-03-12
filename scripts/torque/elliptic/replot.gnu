set terminal pdfcairo enhanced
set output "mx_over_t.pdf"
set xlabel "t [s]"
set ylabel "<m>"
p "m.dat" u 1:2 w lp t "<mx>", "" u 1:3 w lp t "<my>", "" u 1:4 w lp t "<mz>"

while (1){
    replot
    pause .1
}
