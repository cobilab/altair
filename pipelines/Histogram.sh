#!/bin/bash
#
MIN=29200;
MAX=30000;
#
lzma -d -k -c SARS-CoV-2.fa.lzma \
| awk '$0 ~ ">" {print c; c=0;printf substr($0,2,10) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' \
| sed '/^$/d' \
| awk '{print $2}' \
| sort \
| uniq -c \
| awk '{print $1"\t"$2}' > histogram-data.txt;
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
    set xrange [$MIN:$MAX]
    set xtics auto
    set ytics auto
    set grid
    set ylabel "Frequency"
    set xlabel "Sequence length"
    plot "histogram-data.txt" u 2:1 w boxes lc rgb"blue" notitle
EOF
#
