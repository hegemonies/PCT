set terminal png size 800, 700 font 'Times New Roman, 14'

set key top left
set pointsize 2

set ylabel 'Ускорение'
set xlabel 'Количество процессов'

set xzeroaxis lt -1
set yzeroaxis lt -1

set grid xtics lc rgb  '#555555' lw 1 lt 0
set grid ytics lc rgb  '#555555' lw 1 lt 0

set output 'result.png'

plot 'N=100' with linespoints lw 1 lt rgb 'red', \
    'N=400' with linespoints lw 1 lt rgb 'blue', \
    'N=1600' with linespoints lw 1 lt rgb 'green', \
    'N=3000' with linespoints lw 1 lt rgb 'yellow', \
    'N=4000' with linespoints lw 1 lt rgb 'pink', \
    'Линейное' with linespoints lw 1 lt rgb 'black'