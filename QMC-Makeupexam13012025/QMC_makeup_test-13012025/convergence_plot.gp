set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# 1. gnuplot -e "a='B2_rs5'" your_script.gp 

# 2. Construct the filename
fname = a . ".out"

fpng = a . ".png"

set output fpng


# 3. Plot by constructing the command string dynamically
#    %s will be replaced by the content of 'fname'
plot sprintf("< grep elocal %s", fname) u 0:3:4 with errorbars

set output
