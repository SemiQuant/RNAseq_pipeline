#!/bin/bash


still to do
check naming and renaming of files

add in metagenomic search





#RNA pipeline from sputum - host and bacterial
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
  ILLUMINACLIP:"./bin/programs/Trimmomatic-0.36/adapters/TruSeq2-SE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:20

#FastQC
/home/lmbjas002/bin/FastQC/fastqc "${read_file}.trimmed.fq.gz" -o "$out_dir"

# -align to human
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

# 0r STAR
# gen_dir="/researchdata/fhgfs/lmbjas002/Server_Version/RNAseq/HumanGenome/Genome_ref/refs"
# cd $gen_dir
# # Generating genome indexes
# /home/lmbjas002/bin/STAR/bin/Linux_x86_64/STAR \
# --runThreadN $threads \
# --runMode genomeGenerate --genomeDir $gen_dir \
# --genomeFastaFiles $gen_dir/GCF_000001405.34_GRCh38.p8_genomic.fa
# #run alignment
#
# /home/lmbjas002/bin/STAR/bin/Linux_x86_64/STAR \
# --runThreadN $threads \
# --genomeDir $gen_dir \
# --sjdbGTFfile ${gen_dir}/GCF_000001405.34_GRCh38.p8_genomic.gtf \
# --readFilesIn $read_file \
# --readFilesCommand zcat

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
python2.6  /home/lmbjas002/bin/HTSeq/scripts/htseq-count -f bam "accepted_hits.sorted.dedup.bam" $Ref_gtf > "${outdir}/${name}.HTSeq.counts"
# or
#using feature count
/home/lmbjas002/bin/subread/bin/featureCounts -a $Ref_gtf -o "${read_file}.featCount.counts" "accepted_hits.sorted.dedup.bam"

#rename
mv "accepted_hits.sorted.dedup.bam" "${name}_accepted_hits.sorted.dedup.bam"



#align unalined to H37Rv
# -convert unaligned to fq
~/bin/bedtools/bin/bedtools bamtofastq -i "${name}_accepted_hits.sorted.dedup.bam" -fq "${name}.unaligned.fastq"
read_file="${name}.unaligned.fastq"

# -align to h37Rv
/home/lmbjas002/bin/bowtie2-2.2.6/bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$threads" -x "${Script_dir}/refs/H37Rv" -U "$read_file" -S "${read_file}.h37.sam"
  #n-celi number of mismatched bases
  #---score-min set min alignment score

#sort Sam
#sort -k 3,3 -k 4,4n "${read_file}.sam" > "${read_file}.sorted.sam"
#Or to sorted BAM
java -Xmx"$ram"g -jar ~/bin/programs/picard-tools-1.124/picard.jar SortSam \
  INPUT="${read_file}.h37.sam" \
  OUTPUT="${read_file}.h37.bam" \
  SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT

#Index using samtools
/opt/exp_soft/samtools-1.1/samtools index "${read_file}.h37.bam"

#Mark PCR duplicates with PICARD
java -Xmx"$ram"g -jar ~/bin/programs/picard-tools-1.124/picard.jar MarkDuplicates \
        INPUT="${read_file}.h37.bam" \
        OUTPUT="${read_file}.sorted.dedup.h37.bam" \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=TRUE \
        ASSUME_SORTED=TRUE \
        M="${read_file}.dedup.h37.bam.txt"

#Cufflinks
/opt/exp_soft/cufflinks-2.1.1/cufflinks -q -p $threads -o "$out_dir" -m 50 -g "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28.gtf" "${read_file}.sorted.dedup.bam"

#CuffQuant to ref
/opt/exp_soft/cufflinks-2.2.1/cuffquant -q -p $threads -o "$out_dir" "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28.gtf" "${read_file}.sorted.dedup.bam"

#echo "seqname	source	feature	start	end	score	strand	frame	attributes" > "${read_file}.transcripts.gtf"
#grep exon transcripts.gtf >> "${read_file}.exon.transcripts.gtf"

#get sme stats such as number of mapped reads
/opt/exp_soft/samtools-1.1/samtools flagstat "${read_file}.sorted.dedup.bam" > "${read_file}.flagstat.txt"

#get raw counts
python2.6  /home/lmbjas002/bin/HTSeq/scripts/htseq-count -f bam "${read_file}.sorted.dedup.bam" "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28.gtf" > "${read_file}.HTSeq.counts"

#using feature count
/home/lmbjas002/bin/subread/bin/featureCounts -a "${Script_dir}/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28.gtf" -o "${read_file}.featCount.counts" "${read_file}.sorted.dedup.bam"
#might need to be .gff file, not sure

#rename files
mv abundances.cxb "${name}.abundances.cxb"
mv genes.fpkm_tracking "${name}.genes.fpkm_tracking"
mv isoforms.fpkm_tracking "${name}.isoforms.fpkm_tracking"
mv skipped.gtf "${name}.skipped.gtf"
mv transcripts.gtf "${name}.transcripts.gtf"


#delete files
rm *.trimmed.fq.gz
rm *.sam
rm *.fastq.gz.trimmed.fq.gz.bam
rm *.fastq.gz.trimmed.fq.gz.bam.bai

#for further unaligned do metagenomics


