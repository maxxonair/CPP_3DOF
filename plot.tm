#!/usr/bin/gnuplot -persist

set xlabel "Velocity [m/s]"
set ylabel "Thrust [N]"
set grid
set xrange [335:345]
#set yrange [-10:21000] 
#autoscale y
#autoscale y2
set y2tics
set title "TrajSim - trajectory basic "
plot "res.dat" using 1:($26)  axes x1y1, "res.dat" using 1:29  axes x1y2

