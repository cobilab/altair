#!/bin/bash
#
lzma -d -k -c SARS-CoV-2.fa.lzma \
| ./AltaiR filter -a ACGT -min $1 -max $2 > FIL-SARS-CoV-2.fa
#

