set title "static structure factor"
unset xlabel
unset ylabel
set hidden3d
set view 40,50
set ticslevel 0
set contour
set view 45,55
splot[][][0:]'sofk' with lines
