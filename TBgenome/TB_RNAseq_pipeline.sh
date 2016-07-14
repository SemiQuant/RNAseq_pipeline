#!/bin/bash

#input variables
Script_dir=$(dirname "$0")
read_file="${1}"
threads="$2"
ram=$(expr $threads \* 2)
indir=$(dirname "$1")

cd "$indir"
name=$(basename "$1")
out_dir="${indir}/${name%_*}"
mkdir "$out_dir"

#FastQC
#get quality of raw reads
/home/lmbjas002/bin/FastQC/fastqc "$read_file" -o "$out_dir"

#Trim Reads
#sliding window of 4 and keep only if len >= 40
#change adapter if necessary
java -Xmx"$ram"g -jar ~/bin/programs/Trimmomatic-0.32/trimmomatic.jar SE -phred33 \
  "$read_file" \
  "${read_file}.trimmed.fq.gz" \
  ILLUMINACLIP:"~/bin/programs/Trimmomatic-0.32/adapters/TruSeq3-SE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:40

mv "${read_file}.trimmed.fq.gz" "${out_dir}/"
cd "${out_dir}"
read_file="${out_dir}/$name.trimmed.fq.gz"

#Align Reads
/home/lmbjas002/bin/bowtie2-2.2.6/bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$threads" -x "${Script_dir}/refs/H37Rv" -U "$read_file" -S "${read_file}.sam"
#n-celi number of mismatched bases
#---score-min set min alignment score

#sort Sam
#sort -k 3,3 -k 4,4n "${read_file}.sam" > "${read_file}.sorted.sam"
#Or to sorted BAM
java -Xmx"$ram"g -jar ~/bin/programs/picard-tools-1.124/picard.jar SortSam \
  INPUT="${read_file}.sam" \
  OUTPUT="${read_file}.bam" \
  SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT

#Index using samtools
/opt/exp_soft/samtools-1.1/samtools index "${read_file}.bam"

#Mark PCR duplicates with PICARD
java -Xmx"$ram"g -jar ~/bin/programs/picard-tools-1.124/picard.jar MarkDuplicates \
        INPUT="${read_file}.bam" \
        OUTPUT="${read_file}.sorted.dedup.bam" \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE \
        ASSUME_SORTED=TRUE \
        M="${read_file}.dedup.bam.txt"

#Cufflinks
/opt/exp_soft/cufflinks-2.1.1/cufflinks -q -p $threads -o "$out_dir" -m 50 -g "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28_renamed.gtf" "${read_file}.sorted.dedup.bam"

#CuffQuant to ref
/opt/exp_soft/cufflinks-2.2.1/cuffquant -q -p $threads -o "$out_dir" "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28_renamed.gtf" "${read_file}.sorted.dedup.bam"

#echo "seqname	source	feature	start	end	score	strand	frame	attributes" > "${read_file}.transcripts.gtf"
#grep exon transcripts.gtf >> "${read_file}.exon.transcripts.gtf"

#get sme stats such as number of mapped reads
/opt/exp_soft/samtools-1.1/samtools flagstat "${read_file}.sorted.dedup.bam" > "${read_file}.flagstat.txt"

#rename files
mv abundances.cxb "${name}.abundances.cxb"
mv genes.fpkm_tracking "${name}.genes.fpkm_tracking"
mv isoforms.fpkm_tracking "${name}.isoforms.fpkm_tracking"
mv skipped.gtf "${name}.skipped.gtf"
mv transcripts.gtf "${name}.transcripts.gtf"

#zip files
# bgzip "${read_file}.sorted.dedup.bam"

#delete files
rm *.trimmed.fq.gz
rm *.sam
rm *.fastq.gz.trimmed.fq.gz.bam
rm *.fastq.gz.trimmed.fq.gz.bam.bai

