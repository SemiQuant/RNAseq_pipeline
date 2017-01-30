#!/bin/bash
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
/home/lmbjas002/bin/FastQC/fastqc "$read_file" -o "$out_dir"

#Trim Reads
java -Xmx"$ram"g -jar ~/bin/Trimmomatic-0.36/trimmomatic.jar SE -phred33 \
  "$read_file" \
  "${read_file}.trimmed.fq.gz" \
  ILLUMINACLIP:"./bin/programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:20
#fix adapter

mv "${read_file}.trimmed.fq.gz" "$out_dir"
cd "$out_dir"
read_file="${out_dir}/${name}.trimmed.fq.gz"

# Preparing the index for tophat2
cd "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/refs"
# /opt/exp_soft/bowtie2-2.1.0/bowtie2-build "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/ref
# s/GCA_000001405.19_GRCh38.p4_genomic.fa" "/home/lmbjas002/RNA_pipeline/human/pipeline/GenomeRef/refs
# /GCA_000001405.19_GRCh38.p4_genomic"


#Align Reads (Tophat 2)
export PATH=/opt/exp_soft/bowtie2-2.2.6/:$PATH
export PATH=/opt/exp_soft/samtools-1.2/:$PATH

/home/lmbjas002/bin/tophat-2.1.1.Linux_x86_64/tophat2 -o "$out_dir" -p $threads --no-coverage-search  $Ref_bowtie "$read_file"
# /opt/exp_soft/tophat-2.0.10/tophat2 -o "$out_dir" -p $threads --no-coverage-search  $Ref_bowtie "$read_file"

#Index using samtools
/opt/exp_soft/samtools-1.2/samtools index "accepted_hits.bam"

#Mark PCR duplicates with PICARD
java -Xmx"$ram"g -jar ~/bin/picard/picard.jar MarkDuplicates \
        INPUT="accepted_hits.bam" \
        OUTPUT="accepted_hits.sorted.dedup.bam" \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE \
        ASSUME_SORTED=TRUE \
        M="${read_file}.dedup.bam.txt"

#Cufflinks
/opt/exp_soft/cufflinks-2.1.1/cufflinks -q -p $threads -m 50 -b $Ref_name -u -o "$outdir" -g $Ref_gtf "accepted_hits.sorted.dedup.bam"

#CuffQuant to ref
/opt/exp_soft/cufflinks-2.2.1/cuffquant -p $threads $Ref_gtf "accepted_hits.sorted.dedup.bam"

#raw counts
/home/lmbjas002/bin/HTSeq/scripts/htseq-count -f bam "accepted_hits.sorted.dedup.bam" $Ref_gtf > "${outdir}/${name}.HTSeq.counts"

#rename
mv "accepted_hits.sorted.dedup.bam" "${name}_accepted_hits.sorted.dedup.bam"


