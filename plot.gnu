#!/usr/bin/gnuplot -persist

set xlabel "Horizontal Track [km]"
set ylabel "Altitude [km]"
set grid
#set output "traject.png"
#set terminal svg enhanced background rgb 'white'
#set object rectangle from screen 0, 0 to screen 1, 1 behind fillcolor rgb 'white' fillstyle solid nob order
set xrange [350:800]
set yrange [0:80] 
set title "TrajSim - trajectory basic "
plot "res.dat" using ($33 / 1000):($4 / 1000) 

