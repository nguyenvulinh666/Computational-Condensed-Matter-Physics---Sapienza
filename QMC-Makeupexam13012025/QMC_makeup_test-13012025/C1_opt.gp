f(x)=a*(x-beta_opt)
# the value at beta=2 is excluded from the fit because it would require
# a higher order polynomial
fit [4:21]f(x)'estimate_rs.txt'u 1:2:3 via a,beta_opt
set xrange [3:23]
set xlabel'rs'
set ylabel'energy'
set term pdf
set out 'C1_opt.pdf'
plot 'estimate_rs.txt'u 1:2:3 w e,f(x)
set term wxt
rep
