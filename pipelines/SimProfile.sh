#!/bin/bash
#
FIL_DROP=0;
#
W1=1;
W2=1000;
#
cat $1 \
| awk '{ print $1"\t"$2;}' \
| gto_filter -w $W1 -t $2 -d $FIL_DROP > .NCDFiltered.txt
#
cat $1 \
| awk '{ print $1"\t"$2;}' \
| gto_filter -w $W2 -t $2 -d $FIL_DROP > .NCDFiltered2.txt
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "NCDProfile$1.pdf"
    set style line 101 lc rgb '#000000' lt 1 lw 4
    set border 3 front ls 101
    set tics nomirror out scale 0.5
    set format '%g'
    set size ratio 0.2
    set key out horiz center top
    set yrange [$3:$4]
    set xrange [:]
    set xtics auto
    set ytics auto
    set grid
    set ylabel "NCD"
    set xlabel "Time Point"
    set border linewidth 1.5
    set style line 1 lc rgb '#322152' lt 1 lw 1 pt 5 ps 0.4 # --- blue
    set style line 3 lc rgb '#009900' lt 1 lw 1.5 pt 6 ps 0.4 # --- green
    set style line 2 lc rgb '#CC0000' lt 1 lw 1.5 pt 7 ps 0.4 # --- red
    set xtics border in scale 0,0 nomirror rotate by -55  autojustify
    set xtics  norangelimit  font ",8"
    nth(countCol,labelCol,n) = ((int(column(countCol)) % n == 2500) ? stringcolumn(labelCol) : "")
    plot ".NCDFiltered2.txt" using 1:2 title '  w=$W2' with lines ls 2, ".NCDFiltered.txt" using 1:2 title '  w=$W1' with lines ls 1
EOF
#