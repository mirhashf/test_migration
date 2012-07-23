#!/bin/bash

# Separates contigs from a single FASTA file into their own FASTA files.

if [[ $# -lt 1 ]]; then
    "Input FASTA file required."
fi

fasta_file=${1}

csplit -z ${fasta_file} '/^>/' '{*}'

for i in xx*
do
  chr=`head -1 $i | sed 's/>//'`
  echo "Creating FASTA for ${chr}."
  mv $i ${chr}.fa
done
