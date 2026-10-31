rmc
set data, nipl
data x-ray, cua1, exp.nipl, one, -2.0,-2.0,0.0, 2.0,-2.0,0.0, -2.0,2.0,0.0
set mode, 1.0, shift
sele all
set sigma, SIGMA
set cycl, NCYC
LOTSLINE
set disp, 100000
