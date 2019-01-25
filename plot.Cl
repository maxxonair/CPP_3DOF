#!/usr/bin/gnuplot -persist

set xlabel "Velocity [m/s]"
set ylabel "Lift Coefficient [-]"
set grid
set xrange [0:4100]
#set yrange [0:100] 
set title "TrajSim - trajectory basic "
plot "res.dat" using 5:($21) 

