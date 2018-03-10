#!/bin/bash
export PATH=/users/bi/jlimberis/bin/cufflinks-2.2.1.Linux_x86_64:$PATH
export PATH=/users/bi/jlimberis/bin:$PATH
export PATH=/users/bi/jlimberis/bin/bedtools2/bin:$PATH
export PATH=/users/bi/jlimberis/bin/bcftools-1.3.1:$PATH
export PATH=/users/bi/jlimberis/bin/FastQC:$PATH
export PATH=/users/bi/jlimberis/bin/htslib-1.3.2:$PATH
export PATH=/users/bi/jlimberis/bin/STAR-2.5.2b/bin/Linux_x86_64:$PATH
export PATH=/users/bi/jlimberis/bin/subread-1.5.1-Linux-x86_64/bin:$PATH
export PATH=/users/bi/jlimberis/bin/HTSeq-0.6.1/scripts:$PATH
export PATH=/users/bi/jlimberis/.local/bin/:$PATH
export PATH=/users/bi/jlimberis/bin/bowtie2-2.3.0:$PATH
export PATH=/users/bi/jlimberis/bin/samtools-1.3.1:$PATH

TRIM=~/bin/trimmomatic.jar
adapterSE=~/bin/Trimmomatic/adapters/universal.fa
adapterPE=~/bin/Trimmomatic/adapters/TruSeq2-PE.fa
# cut_adapt_seq="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TTACTATGCCGCTGGTGGCTCTAGATGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGA"
PICARD=~/bin/picard.jar
GATK=~/bin/GenomeAnalysisTK.jar
Script_dir=$(dirname "$0")

#check if programs installed
command -v cufflinks >/dev/null 2>&1 || { echo >&2 "I require cufflinks but it's not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo >&2 "I require bcftools but it's not installed. Aborting."; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v STAR >/dev/null 2>&1 || { echo >&2 "I require STAR but it's not installed. Aborting."; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo >&2 "I require htseq but it's not installed. Aborting."; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "I require python2.* but it's not installed. Aborting."; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo >&2 "I require featureCounts but it's not installed. Aborting."; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require bowtie2 but it's not installed. Aborting."; exit 1; }
if [ cut_adapt == "Y"]; then command -v cutadapt >/dev/null 2>&1 || { echo >&2 "I require cutadapt but it's not installed. Aborting."; exit 1; }; fi

if [ ! -f "$TRIM" ]; then echo "$TRIM not found!"; exit 1; fi
if [ ! -f "$PICARD" ]; then echo "$PICARD not found!"; exit 1; fi

if [ $# == 0 ]
  then
    echo -e 'Usage: ./RNA_processes.sh "1 - Input_paramaters.txt" "2 - threads" "3 - ram" \n
    "4 - trim and clip adapt (Y|N)" "5 - is miRNA (Y|N)" "6 - stranded library (yes|no|reverse)" \n
    "7 - cufflinks run in addition to HTseqCount (Y|N)" 8 - "SubRead FeatCount run in addition to HTseqCount (Y|N)" "9 - run qualimap (Y|N)" \n
    "10 - SNP calling from seq data? (Y|N)" \n
    "11 - genome1" "12 - GTF for genome 1" "13 - Type 1 (E=eukaryotic, B=bacterial)" \n
    "14 - genome2" "15 - GTF for genome 2" "16 - Type 2 (E=eukaryotic, B=bacterial)" \n
    if indexing a genome for the first time, this will require >30GB ram for a human genome\n'

    echo -e "Input_paramaters.txt should be a comma seperated list conatining the following:\n
    the directory where the read file(s) are, the output name, the output directory, the fastq file, \n
    the second fastq file if PE reads - else leave blank (i.e file1,) \n"
    exit 1
fi
#/users/bi/jlimberis/CASS_RNAseq,C100,/users/bi/jlimberis/RNAseqData,C100_GTAGAG_HS374-375-376-merged_R1_001.fastq.gz,,/users/bi/jlimberis/testing/Homo_sapiens.GRCh38.dna.primary_assembly.fa,/users/bi/jlimberis/testing/GCF_000195955.2_ASM19595v2_genomic.fna,E,B,/users/bi/jlimberis/testing/Homo_sapiens.GRCh38.87.gtf,/users/bi/jlimberis/testing/GCF_000195955.2_ASM19595v2_genomic.gff

####################
##Define fucntions##
####################

get_reference () {
  mkdir "${Script_dir}/references" #wont overwrite so its ok
  if [[ ! -e $1 ]]
  then
    if [[ ! -e "${Script_dir}/references/$(basename $1).fasta" ]]
    then
        echo "Downloading reference genome $(basename $1)"
        curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$(basename $1)&rettype=fasta" > "${Script_dir}/references/$(basename $1).fasta"
        export ${2}="${Script_dir}/references/$(basename $1).fasta"
    else
        export ${2}="${Script_dir}/$(basename $1).fasta"
        echo "Found reference genome file for $(basename $1)"
    fi
  else
    echo "Found reference genome file for $(basename $1)"
  fi

  if [[ ! -e $3 ]]
  then
    if [[ ! -e "${Script_dir}/references/$(basename $3).gtf" ]]
    then
        curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$(basename $1)&rettype=gtf" > "${Script_dir}/references/$(basename $3).gtf"
        export ${4}="${Script_dir}/references/$(basename $3).gtf"
    else
        export ${4}="${Script_dir}/references/$(basename $3).gtf"
        echo "Found annotations for genome file for $(basename $3)"
    fi
  else
    echo "Found reference genome file for $(basename $1)"
  fi
}

qc_trim_SE () {
  #FastQC pre
  fastqc -t $4 "$1" -o "$2"

  if [[ $trim == "Y" ]]
  then
      if [[ -e "${1/.f*/.trimmed.fq.gz}" ]]
      then
          echo "Found ${1/f*/forward.fq.gz}"
      else
          #Trim Reads
          echo "trimming started $1"
          # java -Xmx"${3}"g -jar ~/bin/trimmomatic.jar SE -phred33 \
          java -jar "$TRIM" SE -phred33 \
            -threads $4 \
            "$1" \
            "${1/.f*/.trimmed.fq.gz}" \
            ILLUMINACLIP:"$adapterSE":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$trim_min

          #FastQC post
          fastqc -t $4 "${1/.f*/.trimmed.fq.gz}" -o "$2"
      fi
  else
      #cheecky workaround after i had alreay written for trmming, this copys the file unnecessarily, will make it better sometime
      cp "$1" "${1/.f*/.trimmed.fq.gz}"

  fi

  echo "trimming completed"


    # if [[ $trim != "Y" ]] && [[ ${#adapterSE} -gt 0 || ${#adapterPE} -gt 0 ]]
    # then
    #
    # fi
}

qc_trim_PE () {
  #FastQC pre
  fastqc -t $5 "$1" -o "$3"
  fastqc -t $5 "$2" -o "$3"

  #Trim Reads
  echo "trimming started $1 $2"
  if [[ $trim == "Y" ]]
  then
      if [[ -e "${1/f*/forward.fq.gz}" ]]
      then
        echo "Found ${1/f*/forward.fq.gz}"
      else
          # java -Xmx"${4}"g -jar ~/bin/trimmomatic.jar PE -phred33 \
          java -jar "$TRIM" PE -phred33 \
            -threads $5 \
            "$1" "$2" \
            "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" \
        		"${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" \
            ILLUMINACLIP:"$adapterPE":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$trim_min

          #FastQC post
          fastqc -t $5 "${1/f*/forward_paired.fq.gz}" -o "$3"
          fastqc -t $5 "${2/f*/_reverse_paired.fq.gz}" -o "$3"

          #as we also want unpaired reads so..
          cat "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" > "${1/f*/forward.fq.gz}"
        	cat "${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" > "${2/f*/reverse.fq.gz}"
      fi
  else
      if [[ -e "${1/f*/forward.fq.gz}" ]]
      then
        echo "Found ${1/f*/forward.fq.gz}"
      else
          #just clip adapeters
          # java -jar "$TRIM" PE -phred33 \
            # -threads $5 \
            # "$1" "$2" \
            # "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" \
        		# "${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" \
            # ILLUMINACLIP:"$adapterPE":2:30:10

          # cat "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" > "${1/f*/forward.fq.gz}"
        	# cat "${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" > "${2/f*/reverse.fq.gz}"
        	cp "$1" "${1/f*/forward.fq.gz}"
        	cp "$2" "${2/f*/reverse.fq.gz}"
      fi
  fi

	echo "trimming completed"


}

BOWTIE_index () {
  #check if indexed alread
  if [ ! -e "${1}.1.bt2" ] #${1/.f*/.1.bt2}
  then
      if [[ $1 == *.gz ]]
      then
          gunzip $1
          bowtie2-build --threads $2 ${1/.gz/} ${1/.gz/} #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      else
          bowtie2-build --threads $2 $1 $1 #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      fi

  else
      echo "Found bowtie2 index for $1?"
  fi
}

STAR_index () {
  #check if indexed alread
  if [[ ! -e "$(dirname $2)/chrLength.txt" ]]
  then
      if [[ $2 == *.gz ]]
      then
        gunzip $2
      fi

      STAR \
      --runThreadN $1 \
      --runMode genomeGenerate --genomeDir $(dirname $2) \
      --genomeFastaFiles "${2/.gz/}" \
      --sjdbGTFfile "$3" \
      --outFileNamePrefix ${2/.f*/}
      # --sjdbOverhang read_length
      #this is the readlength of the RNA data - can get it from fastqs using fastqc or awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' fastq
  else
      echo "Found STAR index for $2?"
  fi
}

BOWTIE_aligner () {
  echo "BOWTIE alignment started $3"
  out_f="${4}/${5}.$(printf $(basename $3) | cut -f 1 -d '.').sam"

  if [[ -e "${out_f/.sam/.bam}" ]]
  then
    echo "Found ${out_f/.sam/.bam}"
  else
    bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -U "$1" -S "$out_f" --un-gz "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
    #$(3 | cut -f 1 -d '.')

    # mv "${4}/un-seqs" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"

    #convert to sorted bam
    # java -Xmx"${6}"g -jar ~/bin/picard*.jar SortSam \
    java -jar "$PICARD" SortSam \
        INPUT="$out_f" \
        OUTPUT="${out_f/.sam/.bam}" \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT

    rm "$out_f"

    # Index using samtools
    # samtools index "${out_f/.sam/.bam}"

    #Mark PCR duplicates with PICARD
    #this is quite slow and not really necessary most of the time
    # java -Xmx"$6"g -jar ~/bin/picard.jar MarkDuplicates \
    #   INPUT="${out_f/.sam/.bam}" \
    #   OUTPUT="${out_f/.bam/.dedup.bam}" \
    #   VALIDATION_STRINGENCY=LENIENT \
    #   REMOVE_DUPLICATES=TRUE \
    #   ASSUME_SORTED=TRUE \
    #   M="${out_f/.bam/.dedup.bam.txt}"
  fi
  echo "BOWTIE alignment completed"
}

BOWTIE_alignerPE () {
  echo "BOWTIE alignment started $3"
  out_f="${4}/${5}.$(printf $(basename $3) | cut -f 1 -d '.').sam"

  if [[ -e "${out_f/.sam/.bam}" ]]
  then
    echo "Found ${out_f/.sam/.bam}"
  else

    bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -1 "$1" -2 "$7" -S "$out_f" --un-gz ${4} --un-conc-gz ${4}
    #$(3 | cut -f 1 -d '.')

    gen=$(basename $3)
    mv "un-conc-mate.1" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
    mv "un-conc-mate.2" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
      # cat "un-seqs" >> xx

    #convert to sorted bam
    # java -Xmx"${6}"g -jar ~/bin/picard*.jar SortSam \
    java -jar "$PICARD" SortSam \
        INPUT="$out_f" \
        OUTPUT="${out_f/.sam/.bam}" \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT

    rm "$out_f"
  fi

  echo "BOWTIE alignment completed"
}

STAR_align () {
  echo "Star alignment started"
   out_f="${4}/${5}.$(printf $(basename $2) | cut -f 1 -d '.').bam"

  if [[ -e "$out_f" ]]
  then
    echo "Found ${out_f}"
  else

    # gtf_file=$(printf $2 | cut -f 1 -d '.')
    gtf_file="$7"
    if [[ $8 == "none" ]]; then
        read2=""
    else
        read2="$8"
    fi

    #use two pass made if intresited in novel jusctions..doubles runtime
    STAR \
        --runThreadN $1 \
        --genomeDir $(dirname $2) \
        --readFilesIn "$3" "$read2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${4}/${5}" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outSAMstrandField intronMotif \
        --sjdbGTFfile $gtf_file
          # --outSAMunmapped

    rm -r "${4}/${5}_STARtmp"
    gen=$(basename $2)
    #ovs this is only needed for PE but doesnt break anything
    mv "${4}/${5}Unmapped.out.mate1" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
    bgzip "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
    mv "${4}/${5}Unmapped.out.mate2" "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"
    bgzip "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"

    mv "${4}/${5}Aligned.sortedByCoord.out.bam" "$out_f"
  fi
  #Index using samtools
  # samtools index "$out_f"

  #Mark PCR duplicates with PICARD
  # java -Xmx"$6"g -jar ~/bin/picard*.jar MarkDuplicates \
  #   INPUT= "$out_f" \
  #   OUTPUT="${out_f/.bam/.dedup.bam}" \
  #   VALIDATION_STRINGENCY=LENIENT \
  #   REMOVE_DUPLICATES=TRUE \
  #   ASSUME_SORTED=TRUE \
  #   M="${out_f/.bam/.dedup.bam.txt}"
  echo "STAR alignment completed"
}

# HISAT_align () {
  # HISAT2 -o "$2" -p $4 --no-coverage-search  $3 "$1"
# }

do_calcs () {
  # gtf_in="$(printf $2 | cut -f 1 -d '.').gtf"
  if [[ $cullfinks == "Y" ]]
  then
    echo "Cufflinks started $4"
    #Cufflinks
    if [[ $read2 == "none" ]]
    then
        cufflinks -q -p $5 -o "$1" -m $7 -g "$4" "$3"
        #-m is average fragment length - ie. for unpaired reads only
    else
        cufflinks -q -p $5 -o "$1" -g "$4" "$3"
    fi

    #CuffQuant to ref
    cuffquant -q -p $5 -o "$1" "$4" "$3"

    #echo "seqname	source	feature	start	end	score	strand	frame	attributes" > "${read_file}.transcripts.gtf"
    #grep exon transcripts.gtf >> "${read_file}.exon.transcripts.gtf"

      #rename files
      mv "${1}/abundances.cxb" "${3/.bam/.abundances.cxb}"
      mv "${1}/genes.fpkm_tracking" "${3/.bam/.genes.fpkm_tracking}"
      mv "${1}/isoforms.fpkm_tracking" "${3/.bam/.isoforms.fpkm_tracking}"
      mv "${1}/skipped.gtf" "${3/.bam/.skipped.gtf}"
      mv "${1}/transcripts.gtf" "${3/.bam/.transcripts.gtf}"
      echo "Cufflinks completed"

  fi

  #get some stats such as number of mapped reads
  #this is outputted by star better but not by bowtie
  samtools flagstat "$3" > "${3/bam/flagstat.txt}"

  # #get raw counts
  # echo "Counts started $4"
  # if [[ $strand == "Y" ]]
  # then
  #   strand2=1
  # elif [[ $strand == "no" ]]
  # then
  #   strand2=0
  # else
  #   strand2=2
  # fi

  if [[ $6 == "B" ]]
  then
      htseq-count --order "pos" --type "gene" -i "Name" --stranded="$strand" -f bam "$3" "$4" > "${3/.bam/.HTSeq.counts}"
      if [[ feat == "Y" ]];then
        featureCounts -t "gene" -g "Name" -O -Q 5 --ignoreDup -T $5 -a "$4" -o "${3/.bam/.featCount.counts}" "$3"
      fi
  # elif [[ $8 == "miRNA" ]]
  # then
      # grep "miRNA" "$4" > "${4/.g*/.miRNA.gtf}"
      # htseq-count --order "pos" --stranded="$strand" -f bam "$3" "${4/.g*/.miRNA.gtf}" > "${3/.bam/.HTSeq.counts}"
      # featureCounts --ignoreDup -T $5 -a "$4" -o "${3/.bam/.featCount.counts}" "$3"
  else
      htseq-count --order "pos" --stranded="$strand" -f bam "$3" "$4" > "${3/.bam/.HTSeq.counts}"
      if [[ feat == "Y" ]];then
        featureCounts --ignoreDup -T $5 -a "$4" -o "${3/.bam/.featCount.counts}" "$3"
      fi
  fi
  echo "Counts completed"

  #can also do qualimap
  export PATH=/users/bi/jlimberis/bin/qualimap_v2.2.1:$PATH
  if [[ $qualimap == "Y" ]]; then
    qualimap rnaseq -bam "$3" -gtf "$4" -outdir  "${3/.bam/_qualimap}"
    # -p "strand-specific-forward"
  fi

}

VaraintCall () {
  #GATK doesnt listen and eats ram so
  jav_ram=$(echo "scale=2; $ram*0.7" | bc)
  export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
    if [ ! -f "$GATK" ]; then
        echo "$GATK not found! Canntor run SNP calling"
    else

    #Add read groups, sort, mark duplicates, and create index
    java -jar $PICARD AddOrReplaceReadGroups \
        I="$2" \
        O="${2}.tmp.snps.bam" \
        SO=coordinate \
        RGID="id" RGLB="library" RGPL="ILLUMINA" RGPU="machine" RGSM="${4}"

        #check if dict exists
        if [ ! -f "${1/.f*/.dict}" ]; then
          java -jar $PICARD CreateSequenceDictionary \
              R="$1"
              O="${1/.f*/.dict}"
        fi

        #check if fai exists
        if [ ! -f "${1}.fai" ]; then
          samtools faidx "$1"
        fi

        /users/bi/jlimberis/ensembl_homosap/Homo_sapiens.GRCh38.dna.fa.fai

        #index bamfile
        samtools index "${2}.tmp.snps.bam"


        # #do this if it compalins
        # java -jar $PICARD ReorderSam \
        #     I="${2}.tmp.snps.bam" \
        #     O="${2}.tmp.snps.reordered.bam" \
        #     R="$1" \
        #     CREATE_INDEX=TRUE
        # rm "${2}.tmp.snps.bam"
        # mv "${2}.tmp.snps.reordered.bam" "${2}.tmp.snps.bam"
        # mv "${2}.tmp.snps.reordered.bai" "${2}.tmp.snps.bam.bai"


    # Split'N'Trim and reassign mapping qualities
        java -jar $GATK \
            -T SplitNCigarReads \
            -R $1 \
            -I "${2}.tmp.snps.bam" \
            -o "${2}.split.bam" \
            -rf ReassignOneMappingQuality \
            -RMQF 255 \
            -RMQT 60 \
            -U ALLOW_N_CIGAR_READS

      rm "${2}.tmp.snps.bam"

      java -jar $PICARD BuildBamIndex \
          I="${2}.split.bam" \
          VALIDATION_STRINGENCY= LENIENT


      #Create a target list of intervals to be realigned with GATK
      java -jar $GATK \
          -T RealignerTargetCreator \
          -R $1 \
          -I "${2}.split.bam" \
          -o "${2}.split.bam.list"
      #-known indels if available.vcf

      #Perform realignment of the target intervals
      java -jar $GATK \
          -T IndelRealigner \
          -R $1 \
          -I "${2}.split.bam" \
          -targetIntervals "${2}.split.bam.list" \
          -o "${2}.tmp2.snps.bam"

        rm "${2}.split.bam"


        # Variant calling
        java -jar $GATK \
            -T HaplotypeCaller \
            -R ${1} \
            -I "${2}.tmp2.snps.bam" \
            -dontUseSoftClippedBases \
            -o "${3}.vcf"

        # rm "${2}.tmp2.snps.bam"
        mv "${2}.tmp2.snps.bam" "${2/.bam/.snps.bam}"


        #Filter - we recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3
        java -jar $GATK \
            -T VariantFiltration \
            -R "${1}" \
            -V "${3}.vcf" \
            -window 35 \
            -cluster 3 \
            -filterName "GATK_recomm" -filter "FS > 30.0 || QD < 2.0" \
            -o "${3}_filtered.vcf"

        #get coverage
        bedtools genomecov -ibam "${2/.bam/.snps.bam}" -bga > "${2/.bam/.bed}"
        bgzip "${2/.bam/.bed}"

    fi

    rm $(ls "${3}.split"*)

    #reset java mem
    jav_ram=$(echo "scale=2; $ram*0.8" | bc)
    export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
}

miRNAaln () {
  #miRNA alignment
  bowtie2 -p $1 --non-deterministic --very-sensitive -x "$2" -U ${3} | samtools view -@ $1 -Sb - > "$4"
  # --very-sensitive = -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
  # or
  # bowtie2 -p "$threads" -D 20 -R 3 -N 1 -L 16 -i S,1,0.50 --non-deterministic -x "$genome1" -U ${read1} | samtools view -@ $threads -Sb - > ${read1/.f*/.bam}
  # or
  # bowtie2 -p "$threads" -N 1 -L 16 --local -x "$genome1" -U ${read1} | samtools view -@ $threads -Sb - > ${read1/.f*/.bam}
  #allow one mismatch for later SNP calling, seed length is 16

  #sort and index
  samtools sort -@ $1 "$4" -o "${4/.bam/.sorted.bam}"
  samtools index "${4/.bam/.sorted.bam}"

  rm "$4"

}

#RNA pipeline = host and bacterial
file_in="$1"
threads="$2"
ram_def=$(expr $threads \* 2)
ram="${3:-$ram_def}"
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
trim="${4:-Y}" #Y|N
trim_min=10
is_mi="${5:-N}"
strand="${6:-reverse}"
cullfinks="${7:-Y}"
feat="${8:-Y}" #subRead feature counts
qualimap="${9:-Y}"
vc="${10:-T}"

genome1="${11:-none}"
G1="${12}"
T1="${13:-E}"
genome2="${14:-none}"
G2="${15}"
T2="${16:-B}"
# cut_adapt="Y"


#look for correct gtf for miRNA else make it
if [[ $is_mi == "Y" ]];
then
  if [[ ! -e ${genome1/.f*/.miRNA"$g_ext"} ]]
    then
      grep "miRNA" $G1 > ${genome1/.f*/.miRNA"$g_ext"}
      G1=${genome1/.f*/.miRNA"$g_ext"}
  fi
  if [[ ! -e ${genome2/.f*/.miRNA"$g_ext"} ]]
    then
    if [[ $genome2 != "none" ]]; then
      grep "miRNA" $G2 > ${genome2/.f*/.miRNA"$g_ext"}
      G2=${genome2/.f*/.miRNA"$g_ext"}
      fi
  fi
fi

while IFS=$',' read -r -a input_vars
do
    read_dir="${input_vars[0]}"
    out_dir="${input_vars[2]:-read_dir}"
    read1="$read_dir/${input_vars[3]}"
    name="${input_vars[1]:-$(basename $read1)}"
    # name="${input_vars[1]:-${input_vars[3]/.f*/}}
    if [[ -z ${input_vars[4]+x} ]]; then read2="$read_dir/${input_vars[4]}"; else read2="none"; fi
    #this isnt working if no comma after so..
    if [[ ${read2} == "$read_dir/" ]]; then read2="none"; fi

    # genome1="${input_vars[5]:-none}"
    # genome2="${input_vars[6]:-none}"
    # T1="${input_vars[7]:-E}"
    # T2="${input_vars[8]:-B}"
    # strand="${input_vars[9]:-reverse}"
    # f_ext="${input_vars[10]:-.fasta}"
    # g_ext="${input_vars[11]:-.gtf}"

    mkdir "${out_dir}/${name}"
    out_dir="${out_dir}/${name}"

    if [[ $genome1 == "none" ]]; then echo "No input genome supplied!"; exit 1; fi

    # G1=${genome1/.f*/$g_ext}
    # G2=${genome2/.f*/$g_ext}

    #set references
    if [[ $genome1 != "none" ]]; then
        get_reference "$genome1" "genome1" "$G1" "G1"; fi
    if [[ $genome2 != "none" ]]; then
        get_reference "$genome2" "genome2" "$G2" "G2"; fi

#####MIRNA - will only be SE so no need for PE options
  if [[ $is_mi == "Y" ]]
  then
    BOWTIE_index $genome1 $threads $G1
    qc_trim_SE "$read1" "$out_dir" $ram $threads
    mv "${read1/.f*/.trimmed.fq.gz}" "$out_dir"
    read1="${out_dir}/$(basename ${read1/.f*/.trimmed.fq.gz})"
    miRNAaln $threads $genome1 $read1 "${out_dir}/${name}.$(basename $genome1).bam"
    bam_file="${out_dir}/${name}.$(printf $(basename $genome1) | cut -f 1 -d '.').bam"
    read_length=$(zcat $read1 | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
    do_calcs $out_dir $genome1 $bam_file $G1 $threads $T1 $read_length "miRNA"

    if [[ $vc = "Y" ]]; then
      VaraintCall "$genome1" "$bam_file" "${out_dir}/${name}" "${name}"
    fi

    if [[ $genome2 != "none" ]]
    then
      BOWTIE_index $genome2 $threads $G2
      mv "${out_dir}/${name}_Unmapped.out.mate1.fastq.gz" "${out_dir}/${name}_${genome1}_Unmapped.out.mate1.fastq.gz"
      gen=$(basename $genome1)
      read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
      miRNAaln $threads $genome2 $read1_unaligned "${out_dir}/${name}.$(basename $genome2).bam"
      bam_file="${out_dir}/${name}.$(printf $(basename $genome2) | cut -f 1 -d '.').bam"
      do_calcs $out_dir $genome2 $bam_file $G2 $threads $T2 $read_length

      if [[ $vc = "Y" ]]; then
        VaraintCall "$genome2" "$bam_file" "${out_dir}/${name}" "${name}"
      fi
    fi

#####MIRNA end
  else
    if [[ $genome1 != "none" ]] && [ $T1 == "E" ]
    then
        STAR_index $threads $genome1 $G1
    elif [ $genome1 != "none" ] && [ $T1 == "B" ]; then
        BOWTIE_index $genome1 $threads $G1
    fi

    if [[ $genome2 != "none" ]]
    then
      if [[ $T2 == "E" ]]; then
        STAR_index $threads $genome2 $G2
      elif [[ $T2 == "B" ]]; then
        BOWTIE_index $genome2 $threads $G2
      else
        echo "Error in reference input 2 - no kingdom selected"
        exit 1
      fi
    fi

    # cd "$indir"

    # cut_adaptpters
    # if [ cut_adapt == "Y"]
    # then
    #   cutadapt -a $cut_adapt_seq -o "${out_dir}/${read1/.f*/.clipped.fq.gz}" "$read1"
    #   read1="${out_dir}/${read1/.f*/.clipped.fq.gz}"
    #   if [[ $read2 == "none" ]]
    #   then
    #     cutadapt -a $cut_adapt_seq -o "${out_dir}/${read2/.f*/.clipped.fq.gz}" "$read2"
    #     read2="${out_dir}/${read2/.f*/.clipped.fq.gz}"
    #   fi
    # fi

    #QC and trim fastq files
    if [[ $read2 == "none" ]]
    then
      qc_trim_SE "$read1" "$out_dir" $ram $threads
      mv "${read1/.f*/.trimmed.fq.gz}" "$out_dir"
      # if [ cut_adapt == "Y"]; then
      #   rm "$read1"
      # fi
      read1="${out_dir}/$(basename ${read1/.f*/.trimmed.fq.gz})"
      if [[ $T1 == "B" ]]
      then
          BOWTIE_aligner "${read1}" "$threads" "$genome1" "$out_dir" "$name" "$ram"
      elif [[ $T1 == "E" ]]
      then
          STAR_align "$threads" "$genome1" "${read1}" "$out_dir" "$name" "$ram" "$G1"
      fi
    else
        qc_trim_PE "$read1" "$read2" "$out_dir" $ram $threads
        mv "${1/f*/forward.fq.gz}" "${2/f*/reverse.fq.gz}" "$out_dir"
        # if [ cut_adapt == "Y"]; then
        #   rm "$read1" "$read2"
        # fi
        read1="${out_dir}/$(basename ${read1/.f*/.trimmed.fq.gz})"
        read2="${out_dir}/$(basename ${read2/.f*/.trimmed.fq.gz})"
        #do PE aligne like above here
        if [[ $T1 == "B" ]]
        then
            BOWTIE_alignerPE "${1/f*/forward.fq.gz}" "$threads" "$genome1" "$out_dir" "$name" "$ram" "${2/f*/reverse.fq.gz}"
        elif [[ $T1 == "E" ]]
        then
            STAR_align "$threads" "$genome1" "${read1/.f*/.trimmed.fq.gz}" "$out_dir" "$name" "$ram" "$G1" "${1/f*/forward.fq.gz}" "${2/f*/reverse.fq.gz}"
        fi
    fi

    bam_file="${out_dir}/${name}.$(printf $(basename $genome1) | cut -f 1 -d '.').bam"
    #this takes the first 2500 reads and calculates the read length
    read_length=$(zcat $read1 | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
    do_calcs $out_dir $genome1 $bam_file $G1 $threads $T1 $read_length

    if [[ $vc = "Y" ]]; then
        VaraintCall "$genome1" "$bam_file" "${out_dir}/${name}" "${name}"
    fi

    if [[ $genome2 != "none" ]]
    then
      if [[ $read2 == "none" ]]
      then
        #convert unaligned to fasta - STAR now has this built in :)
        mv "${out_dir}/${name}_Unmapped.out.mate1.fastq.gz" "${out_dir}/${name}_${genome1}_Unmapped.out.mate1.fastq.gz"
        gen=$(basename $genome1)
        read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
        #what if the first alignement was done with bowtie??
        if [[ $T2 == "B" ]]
        then
          BOWTIE_aligner "$read1_unaligned" "$threads" "$genome2" "$out_dir" "$name" "$ram"
        elif [[ $T2 == "E" ]]
        then
          STAR_align "$threads" "$genome2" "$read1_unaligned" "$out_dir" "$name" "$ram"
        fi
      else
        gen=$(basename $genome1)
        read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
        read2_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate2.fastq.gz"
        if [[ $T2 == "B" ]]
        then
            BOWTIE_aligner "$read1_unaligned" "$threads" "$genome2" "$out_dir" "$name" "$ram" "$read2_unaligned"
        elif [ $T2 == "E" ]
        then
            STAR_align "$threads" "$genome2" "$read1_unaligned" "$out_dir" "$name" "$ram" "$read2_unaligned"
        fi
        bam_file="${out_dir}/${name}.$(printf $(basename $genome2) | cut -f 1 -d '.').bam"
        do_calcs $out_dir $genome2 $bam_file $G2 $threads $T2 $read_length
        if [[ $vc = "Y" ]]; then
            VaraintCall "$genome2" "$bam_file" "${out_dir}/${name}" "${name}"
        fi
      fi
    fi
  fi
#cleanup

#see what shoudl be removed, remember to leave those reads unaligned to genome two, may want to balst them or something

done<"$file_in"




#random things
# perl ~/bin/gtf2bed.pl Homo_sapiens.GRCh38.86.gtf > Homo_sapiens.GRCh38.86.bed
# G1_bed=/users/bi/jlimberis/RNAseqData/ens_gen/Homo_sapiens.GRCh38.86.bed
# cd /users/bi/jlimberis/RNAseqData/ens_gen/trimmed/testing/T006
# read_distribution.py -i T006.miSub.fq.gz.sorted.bam -r $G1_bed

#split reads into small and other RNA
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) <= 30 && length($2) >= 10) {print "@"$0} }' > ${i}.miRNA.fq
#bgzip ${i}.miRNA.fq
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) > 30) {print "@"$0} }' > ${i}.mRNA.fq
#bgzip ${i}.mRNA.fq

















