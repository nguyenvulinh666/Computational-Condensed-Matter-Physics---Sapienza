reset
# 1. Set the file format (pngcairo is better quality than standard png)
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# 2. Set the output filename
set output 'sofk.png'

set title "static structure factor"
unset xlabel
unset ylabel
set hidden3d
set view 40,50
set ticslevel 0
set contour
set view 45,55
splot[][][0:]'sofk' with lines

set output
