#!/usr/bin/gnuplot -persist

set xlabel "Velocity [m/s]"
set ylabel "Altitude [km]"
set grid
#set output "traject.png"
#set terminal svg enhanced background rgb 'white'
#set object rectangle from screen 0, 0 to screen 1, 1 behind fillcolor rgb 'white' fillstyle solid nob order
set xrange [0:5100]
#set yrange [0:60] 
set title "TrajSim - trajectory basic "
plot "res.dat" using 5:($34) w l lw 2
#plot "res2.dat" using 5:($4 / 1000) 

