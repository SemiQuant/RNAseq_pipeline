Script_dir=$(dirname "$0") #get directory where script is
F_dir="$1"
Files_IN="$2"
Out_dir="/researchdata/fhgfs/lmbjas002/NewPipeRes"

#changable parameters
Trim=20
Call=10
Emit=1
threads=12
ram=24
Ref_name="${Script_dir}/references/H37Rv/H37Rv.fasta" #$7

IFS=$'\n'
for sample in $(cat $Files_IN); do
  mkdir "${Out_dir}/${sample}"
  Mdata="${Out_dir}/${sample}/Data"
  Temp="${Out_dir}/${sample}/Temp"
  mkdir "$Temp"
  basic_log="${Mdata}/Run_log_Basic.log"
  START="$(date +%s)" #strat time of script run - end at bottom and diff calc

  Reads=

  #fastQC Sequence_Stats
  # echo "$(date) fastQC Sequence_Stats started" >> "$basic_log"
  # nohup sh "${Script_dir}/FastQC.sh" "$Reads" "$SeqStats"
  # echo "$(date) fastQC Sequence_Stats completed" >> "$basic_log"

  # #Trim sequences and clip adapters
  # java -Xmx"$ram"g -jar ~/bin/Trimmomatic-0.36/trimmomatic.jar SE -phred33 -threads "$threads" \
  #     "$Reads" \
  #     "${Reads}.trimmed.fq.gz" \
  #     ILLUMINACLIP:"~/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:20

  #fastQC Sequence_Stats
  # echo "$(date) fastQC Sequence_Stats started" >> "$basic_log"
  # nohup sh "${Script_dir}/FastQC.sh" "${Reads}.trimmed.fq.gz" "$SeqStats"
  # echo "$(date) fastQC Sequence_Stats completed" >> "$basic_log"

  #bam to fastq
  ~/bin/bedtools/bin/bedtools bamtofastq $Reads "${Reads}.fq.gz"
  bgzip "${Reads}.fq"

  #mapping
  ~/bin/bwa/bwa mem -M -t "$threads" -c 100 -R "@RG\\tID:${Reads}\\tSM:${Reads}\\tPL:Illumina" -M -T 50 "$Ref_name" "${Reads}.fq.gz" > "${Reads}_BWA.sam"

  #convert SAM into sorted BAM via picard
  java -Xmx"${ram}"g -jar ~/bin/picard/picard.jar SortSam \
      INPUT="${Reads}_BWA.sam" \
      OUTPUT="${Reads}_BWA.sorted.bam" \
      SORT_ORDER=coordinate \
      VALIDATION_STRINGENCY=LENIENT

  rm "${Reads}_BWA.sam"

  #Index using samtools
  ~/bin/samtools/samtools index "${Reads}_BWA.sorted.bam"

  #Mark PCR duplicates with PICARD
  java -Xmx"${ram}"g -jar ~/bin/picard/picard.jar MarkDuplicates \
      INPUT="${Reads}_BWA.sorted.bam" \
      OUTPUT="${Reads}_BWA.sorted.dedup.bam" \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=TRUE \
      ASSUME_SORTED=TRUE \
      M="${Reads}_BWA.sorted.dedup.bam.txt"

#   #generate index  with picard
#   java -Xmx"${ram}"g -jar ~/bin/picard/picard.jar BuildBamIndex \
#       I="${Reads}_BWA.sorted.dedup.bam" \
#       VALIDATION_STRINGENCY= LENIENT
#
#   #Create a target list of intervals to be realigned with GATK
#   java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
#       -T RealignerTargetCreator \
#       -R "$Ref_name" \
#       -I "${Reads}_BWA.sorted.dedup.bam" \
#       -o "${Reads}_BWA.sorted.dedup.bam.list"
# 	#-known indels if available.vcf
#
#   #Perform realignment of the target intervals (base quality score recalibrations)
#   java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
#       -T IndelRealigner \
#       -R "$Ref_name" \
#       -I "${Reads}_BWA.sorted.dedup.bam" \
#       -targetIntervals "${Reads}_BWA.sorted.dedup.bam.list" \
#       -o "${Reads}_BWA.sorted.dedup.realigned.bam"

  #sort file with samtools
  ~/bin/samtools/samtools sort -@ $threads -m "${ram}g" -O "bam" -T "working" -o  "${Reads}_BWA.sorted.dedup.sorted.bam" "${Reads}_BWA.sorted.dedup.bam"

  #index the sorted file
  ~/bin/samtools/samtools index "${Reads}_BWA.sorted.dedup.sorted.bam"

  #delete some temp files
  rm *.sam
  rm *dedup.bam

  #CollectMultipleMetrics
  java -Xmx"${ram}"g -jar ~/bin/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT="${Reads}_BWA.sorted.dedup.sorted.bam" \
      OUTPUT="${Reads}"

  #samtools flagstats -get mapping statistics
  ~/bin/samtools/samtools flagstat "${Reads}_BWA.sorted.dedup.sorted.bam"  >> "${Reads}.flagstat.txt"

  ~/bin/bedtools/bin/bedtools genomecov -ibam "${Reads}_BWA.sorted.dedup.sorted.bam" -bga > "${Reads}.genome_coverage.txt"

  depth_cov=$(samtools depth "${Reads}_BWA.sorted.dedup.sorted.bam"  |  awk '{sum+=$3} END { print sum/NR}')
  depth_cov=${depth_cov%.*}
  depth_cov=$((depth_cov*2))
  ~/bin/samtools/samtools mpileup -B -Q 23 -d $((depth_cov*3)) -C 50 -ugf "$Ref_name" "${Reads}_BWA.sorted.dedup.sorted.bam"  2> "pileup.log" | ~/bin/bcftools/bcftools call -vc -Ob --threads "$5" > "${Reads}.raw.bcf"
  #--ploidy "X"
  ~/bin/bcftools/bcftools view "${Reads}.raw.bcf" | ~/bin/bcftools/vcfutils.pl varFilter -d 1 -D $depth_cov > "${Reads}_raw_variants_Samtools_cons_caller.vcf"
  rm "*raw.bcf"
  rm "pileup.log"

  #call variants with GATK
  module load java/jdk-1.7
  java -Xmx"${ram}"g -jar ~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar \
      -nt "${threads}" \
      -T UnifiedGenotyper \
      -R "${Ref_name}" \
      -I "${Reads}_BWA.sorted.dedup.sorted.bam" \
      -ploidy 1 \
      -glm BOTH \
      -stand\_call\_conf "${call}" \
      -stand\_emit\_conf "${emit}" \
      -o "${Reads}_raw_variants_GATK_unified_genotyper.vcf"

cd "$Temp"

done