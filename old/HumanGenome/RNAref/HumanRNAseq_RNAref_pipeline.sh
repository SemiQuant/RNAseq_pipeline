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
#/home/lmbjas002/bin/FastQC/fastqc "$read_file" -o "$out_dir"

#Trim Reads
#java -Xmx"$ram"g -jar ~/bin/programs/Trimmomatic-0.32/trimmomatic.jar SE -phred33 \
#  "$read_file" \
#  "${read_file}.trimmed.fq.gz" \
#  ILLUMINACLIP:"/home/lmbjas002/bin/programs/Trimmomatic-0.32/adapters/TruSeq3-PE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:10

mv "${read_file}.trimmed.fq.gz" "$out_dir"
cd "$out_dir"
#read_file="${out_dir}/${name}.trimmed.fq.gz"

#Align Reads (Tophat 2)
/home/lmbjas002/bin/tophat-2.1.0.Linux_x86_64/tophat2 -o "$out_dir" -p $threads --no-coverage-search "$Script_dir/refs/GCF_000001405.30_GRCh38.p4_rna" "$read_file"
#make sure its using bowtie2

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
/opt/exp_soft/cufflinks-2.1.1/cufflinks -q -p $threads -m 50 -b "${Script_dir}/refs/GCF_000001405.30_GRCh38.p4_rna.fa" -u -o "$out_dir" -g "${Script_dir}/refs/Homo_sapiens.GRCh38.81.gtf" "accepted_hits.sorted.dedup.bam"

#CuffQuant to ref
/opt/exp_soft/cufflinks-2.2.1/cuffquant -q -p $threads -o "$out_dir" "${Script_dir}/refs/Homo_sapiens.GRCh38.81.gtf" "accepted_hits.sorted.dedup.bam"

#zip files
bgzip accepted_hits.sorted.dedup.bam
bgzip unmapped.bam

#rename files
mv abundances.cxb "${name}.abundances.cxb"
mv accepted_hits.sorted.dedup.bam.gz "${name}.accepted_hits.sorted.dedup.bam.gz"
mv align_summary.txt "${name}.align_summary.txt"
mv deletions.bed "${name}.deletions.bed"
mv genes.fpkm_tracking "${name}.genes.fpkm_tracking"
mv insertions.bed "${name}.insertions.bed"
mv isoforms.fpkm_tracking "${name}.isoforms.fpkm_tracking"
mv junctions.bed "${name}.junctions.bed"
mv prep_reads.info "${name}.prep_reads.info"
mv skipped.gtf "${name}.skipped.gtf"
mv transcripts.gtf "${name}.transcripts.gtf"
mv unmapped.bam.gz "${name}.unmapped.bam.gz"

#delte files
rm *trimmed*
rm accepted_hits.bam
rm accepted_hits.bam.bai

