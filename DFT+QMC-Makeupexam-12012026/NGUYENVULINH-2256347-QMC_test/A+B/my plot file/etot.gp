
reset
# 1. Set the file format (pngcairo is better quality than standard png)
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# 2. Set the output filename
set output 'etot.png'

set title "energy per particle"
set xlabel "block index"
set ylabel "etot"
unset key
plot'< grep etot worm.out'u 0:1 with lines,''using 0:3:4 with errorbars

set output
