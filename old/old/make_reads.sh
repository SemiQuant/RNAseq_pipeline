#!/bin/bash
length=20
cd "/Users/jdlim/Desktop/RNA/TB"
ref=$(<H37Rv.fasta.oneline.txt)
IFS=$';\n'
while read -r gene start end
do
  start="$(echo "$start" | tr -d '[[:space:]]')"
  end="$(echo "$end" | tr -d '[[:space:]]')"
  while [ $start -lt $(( $end - $length )) ]
  do
    echo ">${gene}_${length}.${start}" >> ${gene}_${length}.fasta
    echo ${ref:$(( start - 1 )):$length} >> ${gene}_${length}.fasta
#    cut -c $start-$(( $start + $length )) H37Rv.fasta.oneline.txt >> ${gene}_${length}.fasta
    let start=$start+1
  done
done<RefSeq.cds_H37Rv.csv
