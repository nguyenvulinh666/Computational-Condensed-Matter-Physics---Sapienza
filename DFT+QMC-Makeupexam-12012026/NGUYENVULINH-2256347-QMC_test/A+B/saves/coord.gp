unset key
set size square
set title "trace of last stored configuration"
set xlabel "x (Angstrom)"
set ylabel "y (Angstrom)"
plot'restart.coord'using 2:3 with points pointtype 7 pointsize 1.0
