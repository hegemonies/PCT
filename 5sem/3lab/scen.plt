set terminal png size 500, 350 font 'Verdana, 10'
set output '28000.png'
plot 'res_28000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple', \
 'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'
set output '45000.png'
plot 'res_45000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple', \
 'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'
