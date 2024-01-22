#!/bin/bash
#
WINDOW="$2";
OUTPUT_FILE="Plot.pdf"
#
./AltaiR frequency -v -a ACTG < $1 | awk '{ print $1"\t"$2"\t"$3"\t"$4;}' > tmp.data;
./AltaiR average -v -c 1 -p -w $WINDOW < tmp.data > freq-data1.csv
./AltaiR average -v -c 2 -p -w $WINDOW < tmp.data > freq-data2.csv
./AltaiR average -v -c 3 -p -w $WINDOW < tmp.data > freq-data3.csv
./AltaiR average -v -c 4 -p -w $WINDOW < tmp.data > freq-data4.csv
#
./AltaiR nc -v --dna --threads 4 -m 2:1:0:1:0:0.8/0:0:0 -m 3:1:0:0:0:0.9/0:0:0 $1 \
| awk '{ print $1"\t"$2;}' | ./AltaiR average -v -c 2 -p -w $WINDOW > nc-data.csv
#
cat $1 \
| awk '$0 ~ ">" {print c; c=0;printf substr($0,2,10) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' \
| sed '/^$/d' \
| awk '{print $2}' | ./AltaiR average -v -c 1 -p -w $WINDOW > size-data.csv
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Arial,16' size 14in, 11in # Increase the plot size
    set output "$OUTPUT_FILE"
    # Significantly increase spacing between subplots
    set multiplot layout 6, 1 spacing 0, 2
    # Adjust margins to accommodate the larger plot size
    set lmargin 15
    set rmargin 15
    set tmargin 2
    set bmargin 2
    unset key
    # Increase the number of decimal points
    set ytics autofreq
    #
    set ylabel "Length"
    set style line 1 lc rgb '#1c626e' lt 1 lw 3
    plot "size-data.csv" using 0:2 with lines title 'Length' ls 1
    #
    set format y "%.4f"
    set ylabel "NC"
    set style line 1 lc rgb '#ab0d5f' lt 1 lw 3
    plot "nc-data.csv" using 0:2 with lines title 'NC' ls 1
    # Plot for Symbol A
    set ylabel "Frequency A"
    set style line 1 lc rgb '#1E90FF' lt 1 lw 3
    plot "freq-data1.csv" using 0:2 with lines title 'Symbol A' ls 1
    # Plot for Symbol C
    set ylabel "Frequency C"
    set style line 1 lc rgb '#32CD32' lt 1 lw 3
    plot "freq-data2.csv" using 0:2 with lines title 'Symbol C' ls 1
    # Plot for Symbol G
    set ylabel "Frequency G"
    set style line 1 lc rgb '#FFA500' lt 1 lw 3
    plot "freq-data3.csv" using 0:2 with lines title 'Symbol G' ls 1
    # Plot for Symbol T
    set ylabel "Frequency T"
    set style line 1 lc rgb '#FF4500' lt 1 lw 3
    set xlabel "Time"
    plot "freq-data4.csv" using 0:2 with lines ls 1
    unset multiplot
EOF
