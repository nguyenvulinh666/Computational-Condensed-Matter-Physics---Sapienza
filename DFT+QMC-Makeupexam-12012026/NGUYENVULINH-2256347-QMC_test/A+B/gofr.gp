reset
# 1. Set the file format (pngcairo is better quality than standard png)
set terminal pdfcairo size 8in,6in enhanced font 'Verdana,10'

# 2. Set the output filename
set output 'gofr.pdf'

unset key
set title "pair distribution function"
set xlabel "r (Angstrom)"
set ylabel "g(r)"
plot'gofr'using 1:2 with lines,''using 1:3:4 with errorbars lw 0.5

set output
