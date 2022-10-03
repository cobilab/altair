#!/bin/bash
#
lzma -d -k -c SARS-CoV-2.fa.lzma \
| ./AltaiR filter -a ACGT -min $1 -max $2 -i "|2019|" -i "|2020|" -i "|2021|" -i "|2022|" | ./AltaiR filter -p "|2019" -p "|2020" -p "|2021" -p "|2022" -p "2019|" -p "2020|" -p "2021|" -p "2022|"  > FIL-SARS-CoV-2.fa
#

