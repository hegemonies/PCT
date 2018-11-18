set terminal png size 500, 350 font 'Verdana, 10'

set output 'result_800.png'
plot 'total_800.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'result_3200.png'
plot 'total_3200.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'result_16000.png'
plot 'total_16000.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'