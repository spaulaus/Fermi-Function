reset
load '~/.gnuplot'
set border lw 3

unset key

file1 = 'fermi.dat'
#file1 = 'fermiComp-royer-fsCalc.dat'
file2 = 'fermi-rose.dat'

set ylabel "Percent Difference (%)"
set xlabel "Kinetic Energy (keV)" 

plot '<paste fermi.dat fermi-rose.dat' u 2:(($3/$6-1)*100) pt 5 lc rgb prim3 ps 1.5

set terminal postscript eps enhanced color solid "Helvetica,18"
set output 'fermiComp-rose.eps'
replot