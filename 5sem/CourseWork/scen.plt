set terminal png size 500, 350 font 'Times New Roman, 10'
set key top left

set output 'result.png'
plot 'result_100.txt' using 1:2 with linespoints lw 1 lt rgb 'red', \
    'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'
