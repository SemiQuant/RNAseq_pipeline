#!/bin/bash

#using 34gb ram to make H37Rv reference

genomeDir="/home/lmbjas002/SNPrnaSeq/ref"
threads=12
ram=24
genome_name="H37Rv.fasta"
readFilesIn1="/home/lmbjas002/SNPrnaSeq/ERR262978_ERP002485.fastq.gz"
# readFilesIn2=""
outdir="/home/lmbjas002/SNPrnaSeq/test_out"
sample_name="test"

# Only once
# #make index for STAR
# cd "$genomeDir"
# STAR --runMode genomeGenerate --genomeDir "$genomeDir" --genomeFastaFiles "$genome_name" --runThreadN $threads --genomeSAindexNbases 5
# #try this for human
# # STAR --runMode genomeGenerate --genomeDir "$genomeDir" --genomeFastaFiles "$genome_name" --runThreadN $threads --genomeSAindexNbases 12
#
# #make indek for GATK
# java -Xmx"${ram}"g -jar /home/lmbjas002/bin/picard/picard.jar CreateSequenceDictionary \
# R="$genome_name" \
# O="${genome_name/.fasta/}.dict"
#
# samtools faidx "$genome_name"


# Alignment
cd "$outdir"
runDir="1pass"
mkdir $runDir
cd $runDir
STAR --readFilesCommand zcat --genomeDir "$genomeDir" --readFilesIn "$readFilesIn1" --runThreadN $threads  #"$readFilesIn2"
# --genomeLoad LoadAndExit

# 2-pass STAR
genomeDir_2="${genomeDir}/hg19_2pass"
mkdir "$genomeDir_2"
cd "$genomeDir"
STAR --runMode genomeGenerate --genomeDir "$genomeDir_2" --genomeFastaFiles "$genome_name" \
    --sjdbFileChrStartEnd "${outdir}/1pass/SJ.out.tab" --sjdbOverhang 75 --runThreadN $threads

# final alignments
cd "$outdir"
runDir="2pass"
mkdir $runDir
cd $runDir
STAR --readFilesCommand zcat --genomeDir "$genomeDir_2" --readFilesIn "$readFilesIn1" --runThreadN $threads

#cleanup
rm -r "$genomeDir_2" "${outdir}/1pass"

#Add read groups, sort, mark duplicates, and create index
cd "${outdir}/2pass"
java -Xmx"${ram}"g -jar /home/lmbjas002/bin/picard/picard.jar AddOrReplaceReadGroups \
    I="${outdir}/2pass/Aligned.out.sam" \
    O="${sample_name}.bam" \
    SO=coordinate \
    RGID="id" RGLB="library" RGPL="platform" RGPU="machine" RGSM=${sample_name}

rm "${outdir}/2pass/Aligned.out.sam"

java -Xmx"${ram}"g -jar /home/lmbjas002/bin/picard/picard.jar MarkDuplicates \
    I="${sample_name}.bam" \
    O="${sample_name}.dedupped.bam" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M="${sample}.output.metrics"

rm "${sample_name}.bam"

module load java/jdk-1.7
# Split'N'Trim and reassign mapping qualities
java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R ${genomeDir}/${genome_name} \
    -I ${sample_name}.dedupped.bam \
    -o ${sample_name}.dedupped.split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

rm ${sample_name}.dedupped.bam

# Indel Realignment (optional)
#?/

# Base Recalibration (optional)
#?/

# Variant calling
java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R ${genomeDir}/${genome_name} \
    -I ${sample_name}.dedupped.split.bam \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
    -stand_emit_conf 20.0 \
    -o ${sample_name}.vcf

bgzip ${sample_name}.dedupped.split.bam

#Filter
# e recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3
java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R ${genomeDir}/${genome_name} \
    -V ${sample_name}.vcf \
    -window 35 \
    -cluster 3 \
    -filterName FS \
    -filter "FS > 30.0" \
    -filterName QD \
    -filter "QD < 2.0" \
    -o ${sample_name}_filtered.vcf

module load java/jdk-1.8













