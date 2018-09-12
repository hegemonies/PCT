set term png\
    size 1200, 900\
    enhanced\
    font 'Arial, 14'

set output 'result.png'

set xzeroaxis lt -1
set yzeroaxis lt -1

set grid xtics lc rgb  '#555555' lw 1 lt 0
set grid ytics lc rgb  '#555555' lw 1 lt 0

set xtics axis
set ytics axis

set key left top

set multiplot \
	title "Runge's function" \
	layout 2, 2 \
	scale 1.05 \
	offset 0.01, 0 \

set title 'Threashold 10'
plot 'result/threashold10_N1000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'N 1000',\
    'result/threashold10_N10000.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'N 10000',\
    'result/threashold10_N100000.txt' using 1:2 with linespoints lw 1 lt rgb 'green' title 'N 100000',\
    'result/threashold10_N1000000.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'N 1000000',\
    'result/threashold10_N10000000.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'N 10000000',\
    'result/threashold10_N100000000.txt' using 1:2 with linespoints lw 1 lt rgb 'brown' title 'N 100000000',\
    'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'orange' title 'linear'

set title 'Threashold 100'
plot 'result/threashold100_N1000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'N 1000',\
    'result/threashold100_N10000.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'N 10000',\
    'result/threashold100_N100000.txt' using 1:2 with linespoints lw 1 lt rgb 'green' title 'N 100000',\
    'result/threashold100_N1000000.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'N 1000000',\
    'result/threashold100_N10000000.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'N 10000000',\
    'result/threashold100_N100000000.txt' using 1:2 with linespoints lw 1 lt rgb 'brown' title 'N 100000000',\
    'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'orange' title 'linear'

set title 'Threashold 1000'
plot 'result/threashold1000_N1000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'N 1000',\
    'result/threashold1000_N10000.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'N 10000',\
    'result/threashold1000_N100000.txt' using 1:2 with linespoints lw 1 lt rgb 'green' title 'N 100000',\
    'result/threashold1000_N1000000.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'N 1000000',\
    'result/threashold1000_N10000000.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'N 10000000',\
    'result/threashold1000_N100000000.txt' using 1:2 with linespoints lw 1 lt rgb 'brown' title 'N 100000000',\
    'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'orange' title 'linear'

set title 'Threashold 10000'
plot 'result/threashold10000_N1000.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'N 1000',\
    'result/threashold10000_N10000.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'N 10000',\
    'result/threashold10000_N100000.txt' using 1:2 with linespoints lw 1 lt rgb 'green' title 'N 100000',\
    'result/threashold10000_N1000000.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'N 1000000',\
    'result/threashold10000_N10000000.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'N 10000000',\
    'result/threashold10000_N100000000.txt' using 1:2 with linespoints lw 1 lt rgb 'brown' title 'N 100000000',\
    'linear_sheet.txt' using 1:2 with linespoints lw 1 lt rgb 'orange' title 'linear'