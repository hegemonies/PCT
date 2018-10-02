set terminal png size 500, 350 font 'Verdana, 10'
set output 'res_runge.png'
plot 'res_runge.txt' using 1:2 with linespoints lw 1 lt rgb 'purple', \
 'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'
set output 'res_monte.png'
plot 'res_monte.txt' using 1:2 with linespoints lw 1 lt rgb 'purple', \
 'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'