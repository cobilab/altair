#!/bin/bash
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "Histogram.pdf"
    set style line 101 lc rgb '#000000' lt 1 lw 2 
    set border 3 front ls 101
    set tics nomirror out scale 0.75
    set key fixed right top vertical Right noreverse noenhanced autotitle nobox
    set style histogram clustered gap 1 title textcolor lt -1
    set style data histograms
    set xtics border in scale 0,0 nomirror #rotate by -60  autojustify
    set yrange [:]
    set xrange [$1:$2]
    set xtics auto
    set ytics auto
    set style line 4 lc rgb '#CC0000' lt 2 dashtype '---' lw 4 pt 5 ps 0.4 # --- red
    set grid
    set ylabel "Frequency"
    set xlabel "Sequence length"
    set arrow from 29903, graph 0 to 29903, graph 1 nohead ls 4
    plot "histogram-data.txt" u 2:1 w boxes lc rgb"blue" notitle
EOF

