reset
# 1. Set the file format (pngcairo is better quality than standard png)
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# 2. Set the output filename
set output 'np.png'
set title "number of particles"
set xlabel "block index"
set ylabel "N"
unset key
plot'< grep np worm.out'using 0:1 with lines,''using 0:3:4 with errorbars

set output
