#!/bin/bash
cd "/Users/jdlim/Desktop/RNA/TB" 
length=20
counter=1
while [ $counter -lt $(( 1524 - $length )) ]
do
  echo -e ">Rv0001_${length}.${counter}" >> Rv0001_${length}.fasta
  cut -c $counter-$(( $counter + $length )) H37Rv.fasta.oneline.txt >> Rv0001_${length}.fasta
  let counter=counter+1
done
