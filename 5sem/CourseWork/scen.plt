set terminal png size 500, 350 font 'Verdana, 10'
set xlabel "Процессы"
set ylabel "Ускорение"

set output 'result_100.png'
plot 'result_100.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'result_400.png'
plot 'result_400.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'result_1600.png'
plot 'result_1600.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'result_6400.png'
plot 'result_6400.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'