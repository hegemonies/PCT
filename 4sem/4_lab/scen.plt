set terminal png size 1000, 450 font 'Verdana, 10'
set output 'total.png' 
set title "N-body"
set ylabel 'Num threads'
set xlabel 'Speedup'
set xzeroaxis lt -1
set yzeroaxis lt -1
set key left top

set grid xtics lc rgb  '#555555' lw 1 lt 0
set grid ytics lc rgb  '#555555' lw 1 lt 0

set xtics axis
set ytics axis
plot 'var_1.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'var 1', \
 'var_2.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'var 2', \
 'var_3.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'var 3', \
 'var_4.txt' using 1:2 with linespoints lw 1 lt rgb 'green' title 'var 4', \
 'var_5.txt' using 1:2 with linespoints lw 1 lt rgb 'yellow' title 'var 5', \
 'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'linear'