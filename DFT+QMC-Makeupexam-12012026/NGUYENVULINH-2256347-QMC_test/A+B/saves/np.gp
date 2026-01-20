set title "number of particles"
set xlabel "block index"
set ylabel "N"
unset key
plot'< grep np worm.out'using 0:1 with lines,''using 0:3:4 with errorbars
