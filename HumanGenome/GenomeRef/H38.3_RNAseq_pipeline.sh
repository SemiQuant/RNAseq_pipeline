#!/bin/bash
Script_dir=$(dirname "$0")
read_file="${1}"
threads="$2"
ram=$(expr $threads \* 2)
indir=$(dirname "$1")

name=$(basename "$1")
out_dir=${indir}/${name%_*}
mkdir "$out_dir"

cd "$indir"

#FastQC
/home/lmbjas002/bin/FastQC/fastqc "$read_file" -o "$out_dir"

#Trim Reads
java -Xmx"$ram"g -jar ~/bin/programs/Trimmomatic-0.32/trimmomatic.jar SE -phred33 \
  "$read_file" \
  "${read_file}.trimmed.fq.gz" \
  ILLUMINACLIP:"./bin/programs/Trimmomatic-0.32/adapters/TruSeq3-PE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:10
#fix adapter

mv "${read_file}.trimmed.fq.gz" "$out_dir"
cd "$out_dir"
read_file="${out_dir}/${name}.trimmed.fq.gz"

#Align Reads (Tophat 2)
/opt/exp_soft/tophat-2.0.10/tophat2 -o "$out_dir" -p $threads --no-coverage-search "$Script_dir/refs/GCA_000001405.19_GRCh38.p4_genomic" "$read_file"

#Index using samtools
/opt/exp_soft/samtools-1.1/samtools index "accepted_hits.bam"

#Mark PCR duplicates with PICARD
java -Xmx"$ram"g -jar ~/bin/programs/picard-tools-1.124/picard.jar MarkDuplicates \
        INPUT="accepted_hits.bam" \
        OUTPUT="accepted_hits.sorted.dedup.bam" \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE \
        ASSUME_SORTED=TRUE \
        M="${read_file}.dedup.bam.txt"

#Cufflinks
/opt/exp_soft/cufflinks-2.1.1/cufflinks -q -p $threads -m 50 -b "${Script_dir}/refs/GCF_000001405.30_GRCh38.p4_rna.fa" -u -o "$outdir" -g "${Script_dir}/refs/Homo_sapiens.GRCh38.81.gtf" "accepted_hits.sorted.dedup.bam"

#CuffQuant to ref
/opt/exp_soft/cufflinks-2.2.1/cuffquant -p $threads "${Script_dir}/refs/Homo_sapiens.GRCh38.81.gtf" "accepted_hits.sorted.dedup.bam"
