

reset
# 1. Set the file format (pngcairo is better quality than standard png)
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# 2. Set the output filename
set output 'coord.png'
unset key
set size square
set title "trace of last stored configuration"
set xlabel "x (Angstrom)"
set ylabel "y (Angstrom)"
plot'restart.coord'using 2:3 with points pointtype 7 pointsize 1.0

set output
