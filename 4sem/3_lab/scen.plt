set terminal png size 500, 350 font 'Verdana, 10'
set output 'total.png'
plot 'runge.txt' using 1:2 with linespoints lw 1 lt rgb 'purple', \
 'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'