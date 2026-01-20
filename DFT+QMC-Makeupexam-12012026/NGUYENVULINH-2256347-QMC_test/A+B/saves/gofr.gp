unset key
set title "pair distribution function"
set xlabel "r (Angstrom)"
set ylabel "g(r)"
plot'gofr'using 1:2 with lines,''using 1:3:4 with errorbars
