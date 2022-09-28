#!/bin/bash
#
# bash Mutated.sh ORIGINAL.fa 71231 -e 0.05 > ORIGINAL_MUTATED.fa
#
gto_fasta_mutate -s $2 -e $3 < $1 \
| gto_fasta_to_seq \
| gto_fasta_from_seq -l 70 -n "Mutated seq" 
#
