#! /software/gnuplot/5.2.0b2/bin/gnuplot

set pm3d map
splot '2DES_T0.txt' u 1:2:3
set terminal png crop enhanced
set output '2DES_T0_real.png'

#set palette defined ( 0 0.05 0.05 0.2, 0.1 0 0 1, 0.25 0.7 0.85 0.9,\
#     0.4 0 0.75 0, 0.5 1 1 0, 0.7 1 0 0, 0.9 0.6 0.6 0.6,\
#     1 0.95 0.95 0.95 )

#set palette defined (0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
#                     0.6667 1 0 0, 1 1 1 1 )
set palette defined (-6 "blue", 0 "white", 3 "red")
set size square
set contour base
set cntrparam level incremental -6, 1, 3
set style line  1 lc rgb "black" lw 1
set style line  2 lc rgb "black" lw 1
set style line  3 lc rgb "black" lw 1
set style line  4 lc rgb "black" lw 1
set style line  5 lc rgb "black" lw 1
set style line  6 lc rgb "black" lw 1
set style line  7 lc rgb "black" lw 1
set style line  8 lc rgb "black" lw 1
set style line  9 lc rgb "black" lw 1
set style line  10 lc rgb "black" lw 1
set style line  11 lc rgb "black" lw 1
set style line  12 lc rgb "black" lw 1
set style line  13 lc rgb "black" lw 1
set style line  14 lc rgb "black" lw 1
set style line  15 lc rgb "black" lw 1
set style increment user
set xrange [33 to 45]
set yrange [33 to 45]
set cbrange[-6 to 3]
set xlabel 'excitation frequency (cm^{-1} x 10^3)'
set ylabel 'detection frequency (cm^{-1} x 10^3)'
#set isosample 250, 250
set pm3d interpolate 0,0
unset key
replot
