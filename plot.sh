# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'rgb_variable.4.png'
set border -1 front linetype -1 linewidth 1.000
set angles degrees
set ticslevel 0
set xtics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 
set title "Test Graph" 
set xlabel "x axis" 
set xlabel  offset character 0, 0, 0 font "" textcolor rgb "red"  norotate
set xrange [ 0.00000 : 200.000 ] noreverse nowriteback
set ylabel "y axis" 
set ylabel  offset character 0, 0, 0 font "" textcolor rgb "green"  rotate by -270
set yrange [ 0.00000 : 200.000 ] noreverse nowriteback
set zlabel "z axis" 
set zlabel  offset character 0, 0, 0 font "" textcolor rgb "blue"  norotate
set zrange [ 0.00000 : 123.000 ] noreverse nowriteback
set lmargin  5
set bmargin  2
set rmargin  5
set cbrange [0:100]
splot 'data.dat' using 1:2:3:4 with points pointsize 1 pt 7 palette	title "Test Graph"
pause -1
