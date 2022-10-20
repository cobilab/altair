#!/bin/bash
#
ITERATIONS=10000;
RATE=0.0001;
#
rm -f ALL.fa;
#
printf "\n" > DIV;
#
gto_genomic_gen_random_dna -n 5000 -s 7 \
| gto_fasta_from_seq -n "Synthetic DNA" > ORIGINAL.fa
#
cp ORIGINAL.fa IN.fa
#
for((x=1 ; x<=$ITERATIONS ; ++x));
  do
  #
  printf "\rRunning iteration $x ...";
  #
  gto_fasta_mutate -s $x -e $RATE < IN.fa \
  | gto_fasta_to_seq \
  | gto_fasta_from_seq -l 70 -n "Synthetic DNA i=$x" > OUT.fa
  cat OUT.fa DIV >> ALL.fa
  cp OUT.fa IN.fa
  done
#
rm -f TMP.fa
#
echo "";
