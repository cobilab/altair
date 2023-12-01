#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file.csv>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="FreqProfilecomp-data.csv.pdf"

# Gnuplot settings for an improved plot with increased spacing
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Arial,10' size 11in, 8.5in # Increase the plot size
    set output "$OUTPUT_FILE"
    
    # Significantly increase spacing between subplots
    set multiplot layout 4, 1 spacing 0, 2

    # Adjust margins to accommodate the larger plot size
    set lmargin 15
    set rmargin 15
    set tmargin 2
    set bmargin 2

    # Increase the number of decimal points
    set format y "%.4f"
    set ytics autofreq

    # Plot for Symbol A
    set ylabel "Freq of A"
    set style line 1 lc rgb '#1E90FF' lt 1 lw 2
    plot "$INPUT_FILE" using 0:1 with lines title 'Symbol A' ls 1

    # Plot for Symbol C
    set ylabel "Freq of C"
    set style line 1 lc rgb '#32CD32' lt 1 lw 2
    plot "$INPUT_FILE" using 0:2 with lines title 'Symbol C' ls 1

    # Plot for Symbol G
    set ylabel "Freq of G"
    set style line 1 lc rgb '#FFA500' lt 1 lw 2
    plot "$INPUT_FILE" using 0:3 with lines title 'Symbol G' ls 1

    # Plot for Symbol T
    set ylabel "Freq of T"
    set style line 1 lc rgb '#FF4500' lt 1 lw 2
    set xlabel "Time"
    plot "$INPUT_FILE" using 0:4 with lines title 'Symbol T' ls 1
    
    unset multiplot
EOF
