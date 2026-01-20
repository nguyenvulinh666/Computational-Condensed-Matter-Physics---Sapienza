set title "energy per particle"
set xlabel "block index"
set ylabel "N"
unset key
plot'< grep etot worm.out'u 0:1 with lines,''using 0:3:4 with errorbars
