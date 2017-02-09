#!/bin/bash

#more than one alignemnt from same sample fastq file?

#








Script_dir=$(dirname "$0")
read_file="${1}"
threads="$2"
ram=$(expr $threads \* 2)
indir=$(dirname "$1")
name=$(basename "$1")
out_dir=${indir}/${name%_*}
mkdir "$out_dir"
Ref_name="${Script_dir}/refs/GCF_000001405.34_GRCh38.p8_genomic.fna"
Ref_gtf="${Script_dir}/refs/GCF_000001405.34_GRCh38.p8_genomic.gff"
Ref_bowtie="Script_dir/refs/GCF_000001405.34_GRCh38.p8_genomic"

cd "$indir"

#FastQC
fastqc "$read_file" -o "$out_dir"

#Trim Reads
java -Xmx"$ram"g -jar Trimmomatic-0.36/trimmomatic.jar SE -phred33 \
  "$read_file" \
  "${read_file}.trimmed.fq.gz" \
  ILLUMINACLIP:"Trimmomatic-0.36/TruSeq3-PE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:20


mv "${read_file}.trimmed.fq.gz" "$out_dir"
cd "$out_dir"
read_file="${out_dir}/${name}.trimmed.fq.gz"

# Preparing the index for tophat2
cd "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/refs"
# /opt/exp_soft/bowtie2-2.1.0/bowtie2-build "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/ref
# s/GCA_000001405.19_GRCh38.p4_genomic.fa" "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/refs
# /GCA_000001405.19_GRCh38.p4_genomic"


#Align Reads (Tophat 2)
tophat2 -o "$out_dir" -p $threads --no-coverage-search  "$Ref_bowtie" "$read_file"

#Index using samtools
samtools index "accepted_hits.bam"

#Mark PCR duplicates with PICARD
java -Xmx"$ram"g -jar picard.jar MarkDuplicates \
        INPUT="accepted_hits.bam" \
        OUTPUT="accepted_hits.sorted.dedup.bam" \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE \
        ASSUME_SORTED=TRUE \
        M="${read_file}.dedup.bam.txt"

#Cufflinks
cufflinks -q -p $threads -m 50 -b $Ref_name -u -o "$outdir" -g $Ref_gtf "accepted_hits.sorted.dedup.bam"

#CuffQuant to ref
cuffquant -p $threads $Ref_gtf "accepted_hits.sorted.dedup.bam"

#raw counts
htseq-count -f bam "accepted_hits.sorted.dedup.bam" $Ref_gtf > "${outdir}/${name}.HTSeq.counts"

#rename
mv "accepted_hits.sorted.dedup.bam" "${name}_accepted_hits.sorted.dedup.bam"

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

#delete files
rm *trimmed*
rm accepted_hits.bam
rm accepted_hits.bam.bai


