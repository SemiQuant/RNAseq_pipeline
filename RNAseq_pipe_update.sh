#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  -dl|--container = download the singularity container to this path and exit
  -t|--threads
  -g1|--genome_reference1 = path to genome reference 1, if only a name is supplied the file will be downloaded from ncbi
  -g2|--genome_reference2 (optional) = path to genome reference 2, if only a name is supplied the file will be downloaded from ncbi
  -g1|--genome_reference1
  -g2|--genome_reference2
  -gtf1|--GTF_reference1
  -gtf2|--GTF_reference2
  -t1|--Type_1 = E for eukaryotic or B forbacterial
  -t2|--Type_2 = E for eukaryotic or B forbacterial
  -r|--ram
  -rd|--read_dir
  -r1|--read1
  -r2|--read2 = the second fastq file if PE reads - else leave blank
  -o|--out_dir
  -n|--name
  -m|--miRNA = Is this miRNA?
  -c|--cufflinks = Run cufflinks?
  -f|--feat_count = Run subread feature count?
  -q|--qualimap = Run qualimap?
  -v|--variant_calling = Perform variant callineg?
  -s|--strand = stranded library (yes|no|reverse)
  -tr|--trim = trim reads?
  -sd|--script_directory
  -fq|--fastQC = run fastqc?
  
  -m|--get_metrics = supply a dir and get metrics for all analyses in that dir, all other options will be ignored if this is non-empyt
  

  Notes
    Think the REST for genomes has changed
    Haven't tested everything since changing format of things
    If it finds a gff make it convert to gff (see code notes at EOF);
      ..No need to do that, because GTF is a tightening of the GFF format. Hence, all GTF files are GFF files, too. 
    
    allow user to pass options to programes, like htseqCount
    lot of mixing of global and local variables, cleanup
    add tmp dir
    
    change setting of flags as "Y" to jsut presence or absence
  "
}

if [ $# == 0 ]
then
    usage
    exit 1
fi


declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -t|--threads)
        threads="$2"
        ;;
        -g1|--genome_reference1)
        g1="$2"
        ;;
        -g2|--genome_reference2)
        g2="$2"
        ;;
        -gtf1|--GTF_reference1)
        gt1="$2"
        ;;
        -gtf2|--GTF_reference2)
        gt2="$2"
        ;;
        -t1|--Type_1) #_(E=eukaryotic_B=bacterial)
        t1="$2"
        ;;
        -t2|--Type_2)
        t2="$2"
        ;;
        -r|--ram)
        ram_in="$2"
        ;;
        -rd|--read_dir)
        read_dir="$2"
        ;;
        -r1|--read1) 
        read1="$2"
        ;;
        -r2|--read2) 
        read2="$2"
        ;;
        -o|--out_dir) 
        out_dir="$2"
        ;;
        -n|--name) 
        name="$2"
        ;;
        -m|--miRNA) #Y or N
        is_mi="Y"
        ;;
        -c|--cufflinks) #Y or N
        cullfinks="Y"
        ;;
        -f|--feat_count) #Y or N
        feat="Y"
        ;;
        -q|--qualimap) #Y or N
        qualimap="Y"
        ;;
        -v|--variant_calling) #Y|N
        vc="Y"
        ;;
        -s|--strand) #stranded library (yes|no|reverse)
        strand="$2"
        ;;
        -tr|--trim) #Y|N
        trim="Y"
        ;;
        -sd|--script_directory)
        Script_dir="$2"
        ;;
        -k|--keep_unpaired)
        keep_unpaired="Y"
        ;;
        -c2|--only_care)
        only_care="Y"
        ;;
        -dl|--container)
        container="$2"
        ;;
        -m|--get_metrics)
        get_metrics="$2"
        ;;
        -fq|--fastQC)
        fastQC="Y"
        ;;
    esac
        shift
    done
}

# supply path to genome files, you can leave blanks and only supply name if you wnat them downloaded, if they cannot be found they will be downloaded

get_reference () {
    mkdir "${Script_dir}/references" 2>/dev/null #wont overwrite so its ok
    if [[ ! -e $1 ]] #check if the file they supplied exists
    then
        local g1_f=$(basename $1)
        local g1_f=${g1_f%.f*}
        if [[ ! -e "${Script_dir}/references/${g1_f}/${g1_f}.fasta" ]] #check if the file they supplied exists in the references folder
        then
            mkdir "${Script_dir}/references/${g1_f}"
            echo "Downloading reference genome ${g1_f}"
            curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${g1_f}&rettype=fasta" >  "${Script_dir}/references/${g1_f}/${g1_f}.fasta"
            export ${2}="${Script_dir}/references/${g1_f}/${g1_f}.fasta"
        else
            export ${2}="${Script_dir}/references/${g1_f}/${g1_f}.fasta"
            echo "Found reference genome file for ${g1_f}"
        fi
    else
        echo "Found reference genome file for $(basename $1)"
    fi
  
    if [[ ! -e $3 ]]
    then
        local gt1_f=$(basename $1)
        local gt1_f=${gt1_f%.f*}
        if [[ ! -e "${Script_dir}/references/${gt1_f}/${gt1_f}.gtf" ]]
        then
            mkdir "${Script_dir}/references/${gt1_f}" 2>/dev/null
            echo "Downloading gtf reference ${gt1_f}"
            curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${gt1_f}&rettype=gtf" >  "${Script_dir}/references/${gt1_f}/${gt1_f}.gtf"
            export ${4}="${Script_dir}/references/${gt1_f}/${gt1_f}.gtf"
        else
            export ${4}="${Script_dir}/references/${gt1_f}/${gt1_f}.gtf"
            echo "Found annotations for genome file for ${gt1_f}"
        fi
    else
        echo "Found reference genome file for $(basename $1)"
    fi
}


STAR_index () {
    if [[ ! -e "$(dirname $2)/chrLength.txt" ]] #check if indexed alread, not robust as require each to be in own folder..but so does star file naming..
    then
        if [[ $2 == *.gz ]]
        then
            gunzip $2
        fi
        if [[ $3 == *.gz ]]
        then
            gunzip $3
        fi
    echo "Star indexing requires about 30GB ram for human genome, so if an error then check the log"
    STAR \
      --runThreadN "$1" \
      --runMode genomeGenerate \
      --genomeDir $(dirname $2) \
      --genomeFastaFiles "${2/.gz/}" \
      --sjdbGTFfile "${3/.gz/}" \
      --outFileNamePrefix ${2/.f*/}
      # --sjdbOverhang read_length
      #this is the readlength of the RNA data - can get it from fastqs using fastqc or awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' fastq
    else
        echo "Found STAR index for $2?"
    fi
}


qc_trim_SE () {
      #FastQC pre
      if [[ ! -e "${2}/${1/.f*/_fastqc.zip}" ]] && [[ ! -z $fastQC ]]
      then
          fastqc -t $4 "$1" -o "$2"
      fi
      
      if [[ ! -z $trim ]]
      then
          if [[ -e "${2}/$(basename ${1/.f*/.trimmed.fq.gz})" ]]
          then
              echo "Found ${1/f*/forward.fq.gz}"
          else
              #Trim Reads
              echo "trimming started $1"
              java -jar "$TRIM" SE -phred33 \
              -threads $4 \
              "$1" \
              "${1/.f*/.trimmed.fq.gz}" \
              ILLUMINACLIP:"$5":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$6
              
              #FastQC post
              if  [[ ! -z $fastQC ]]
              then
                  fastqc -t $4 "${1/.f*/.trimmed.fq.gz}" -o "$2"
              fi
          fi
          mv "${1/.f*/.trimmed.fq.gz}" "2"
          # export read1="${2}/$(basename ${1/.f*/.trimmed.fq.gz})"
          read1="${2}/$(basename ${1/.f*/.trimmed.fq.gz})"
          export read1
      fi
      echo "trimming completed"
}


qc_trim_PE () {
    #FastQC pre
    if [[ ! -e "${3}/${1/.f*/_fastqc.zip}" ]] || [[ ! -e "${3}/${2/.f*/_fastqc.zip}" ]]
    then
        if [[ ! -z $fastQC ]]
        then
            fastqc -t $5 "$1" -o "$3"
            fastqc -t $5 "$2" -o "$3"
        fi
    fi
    
    #Trim Reads
    echo "trimming started $1 $2"
    if [[ ! -z $trim ]]
    then
        # if [[ -e "${1/f*/forward.fq.gz}" ]]
        read1="${3}/$(basename ${1/.f*/_forward.fq.gz})"
        read2="${3}/$(basename ${2/.f*/_reverse.fq.gz})"
        if [[ -e "$read1" ]] && [[ -e "$read2" ]]
        then
            echo "Found ${1/.f*/}"
            export read1
            export read2
        else
            java -jar "$TRIM" PE -phred33 \
              -threads $5 \
              "$1" "$2" \
              "${1/.f*/_forward_paired.fq.gz}" "${1/.f*/_forward_unpaired.fq.gz}" \
              "${2/.f*/_reverse_paired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" \
              ILLUMINACLIP:"$6":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$7
        

            # mv "${1/.f*/_forward.fq.gz}" "${2/.f*/_reverse.fq.gz}" "$3"
            mv "${1/.f*/_forward_paired.fq.gz}" "$read1"
            mv "${2/.f*/_reverse_paired.fq.gz}" "$read2"
            
            if [[ ! -z $keep_unpaired ]]
            then
                #as we also want unpaired reads so..
                # woudl be quicker if you do a cp and then a cat
                # cat "${1/.f*/_forward_paired.fq.gz}" "${1/.f*/_forward_unpaired.fq.gz}" > "$read1"
                # cat "${2/.f*/_reverse_paired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" > "$read2"
                cat "${1/.f*/_forward_unpaired.fq.gz}" >> "$read1"
                cat "${2/.f*/_reverse_unpaired.fq.gz}" >> "$read2"
                
                # just to make sure as sometimes it lets one through if merging the paired and unpaired
                cutadapt --cores=$5 --quality-cutoff 10,10 --minimum-length $7 -o "${read1}.tmp.gz" "$read1"
                mv "${read1}.tmp.gz" "$read1"
                cutadapt --cores=$5  --quality-cutoff 10,10 --minimum-length $7 -o "${read2}.tmp.gz" "$read2"
                mv "${read2}.tmp.gz" "$read2"
            # else
            #     mv "${1/.f*/_forward_paired.fq.gz}" "$read1"
            #     mv "${2/.f*/_reverse_paired.fq.gz}" "$read2"
            fi
            
            mv "${1/.f*/_forward_unpaired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" "${3}"
            
            #FastQC post
            if [[ ! -z $fastQC ]]
            then
                fastqc -t $5 "$read1" -o "$3"
                fastqc -t $5 "$read2" -o "$3"
            fi
        fi
        # read1="${3}/${1/.f*/_forward.fq.gz}"
        export read1
        export read2
    # else
        # if [[ -e "${1/f*/forward.fq.gz}" ]]
        # then
        #     echo "Found ${1/f*/forward.fq.gz}"
    # else
      #just clip adapeters
      # java -jar "$TRIM" PE -phred33 \
      # -threads $5 \
      # "$1" "$2" \
      # "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" \
      # "${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" \
      # ILLUMINACLIP:"$adapterPE":2:30:10
      
      # cat "${1/f*/forward_paired.fq.gz}" "${1/f*/_forward_unpaired.fq.gz}" > "${1/f*/forward.fq.gz}"
      # cat "${2/f*/_reverse_paired.fq.gz}" "${2/f*/_reverse_unpaired.fq.gz}" > "${2/f*/reverse.fq.gz}"
        # fi
        # cp "$1" "${1/.f*/_forward.fq.gz}"
        # cp "$2" "${2/.f*/_reverse.fq.gz}"
    fi
    echo "trimming completed"
}


BOWTIE_index () {
  if [ ! -e "${1/.f*/}.1.bt2" ] #check if indexed alread #${1/.f*/.1.bt2}
  then
      if [[ $1 == *.gz ]]
      then
          gunzip $1
          bowtie2-build --threads $2 ${1/.gz/} ${1/.f*/} #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      else
          bowtie2-build --threads $2 $1 ${1/.f*/} #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      fi
  else
      echo "Found bowtie2 index for $1?"
  fi
}


BOWTIE_alignerSE () {
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
    # java -jar "$PICARD" SortSam \
    #   -INPUT "$out_f" \
    #   -OUTPUT "${out_f/.sam/.bam}" \
    #   -SORT_ORDER coordinate \
    #   -VALIDATION_STRINGENCY LENIENT
    java -jar "$PICARD" SortSam \
      I="$out_f" \
      O="${out_f/.sam/.bam}" \
      SORT_ORDER=queryname \
      VALIDATION_STRINGENCY=LENIENT
      
    rm "$out_f"
    fi
    echo "BOWTIE alignment completed"
}


BOWTIE_alignerPE () {
    echo "BOWTIE alignment started $3"
    # BOWTIE_alignerPE "$read1_unaligned" "$threads" "$g2" "$out_dir" "$name" "$ram" "$read2_unaligned"
    out_f="${4}/${5}.$(printf $(basename $3) | cut -f 1 -d '.').sam"
    
    
    # -U <r> # Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the “standard in” or “stdin” filehandle.

    
    
    if [[ -e "${out_f/.sam/.bam}" ]]
    then
        echo "Found ${out_f/.sam/.bam}"
    else
        # bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -1 "$1" -2 "$7" -S "$out_f" --un-gz ${4} --un-conc-gz ${4}
         bowtie2 \
            --dovetail \
            --local \
            --minins 0 \
            -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 \
            -p "$2" -x ${3/.f*/} -1 "$1" -2 "$7" -S "$out_f" --un-gz ${4} --un-conc-gz ${4}
        
    #$(3 | cut -f 1 -d '.')
    gen=$(basename $3)
    mv "un-conc-mate.1" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
    mv "un-conc-mate.2" "${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"
    export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
    export read2_unaligned="${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"
        
    # cat "un-seqs" >> xx
    
    #convert to sorted bam
    # New syntax
    # java -jar "$PICARD" SortSam \
    #   -INPUT "$out_f" \
    #   -OUTPUT "${out_f/.sam/.bam}" \
    #   -SORT_ORDER queryname \
    #   -VALIDATION_STRINGENCY LENIENT
    java -jar "$PICARD" SortSam \
      I="$out_f" \
      O="${out_f/.sam/.bam}" \
      SORT_ORDER=queryname \
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
        if [[ $(basename $8) == "none" ]]; then
            read2=""
        else
            read2="$8"
        fi
        #use two pass mode if intresited in novel jusctions..doubles runtime
        STAR \
          --runThreadN $1 \
          --genomeDir $(dirname $2) \
          --readFilesIn "$3" "$read2" \
          --readFilesCommand zcat \
          --outFileNamePrefix "${4}/${5}" \
          --outSAMtype BAM Unsorted \
          --outReadsUnmapped Fastx \
          --outSAMstrandField intronMotif \
          --sjdbGTFfile "$7" \
          --quantMode GeneCounts #The counts coincide with those produced by htseq-count with default parameters. 
          # --outSAMunmapped
    
        mv "${5}ReadsPerGene.out.tab" "${4}/${5}_ReadsPerGene.out.tab"
        # rm -r "${4}/${5}_STARtmp"
        gen=$(basename $2)
        #ovs this is only needed for PE but doesnt break anything
        mv "${4}/${5}Unmapped.out.mate1" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
        bgzip "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
        mv "${4}/${5}Unmapped.out.mate2" "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"
        bgzip "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"
        
        export read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
        export read2_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate2.fastq.gz"
        
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


do_calcs () {
  # do_calcs $out_dir $g1 $bam_file $gt1 $threads $t1 $read_length
    # gtf_in="$(printf $2 | cut -f 1 -d '.').gtf"
    
    
    # For counting PE fragments associated with genes, the input bam files need to be sorted by read name 
    # (i.e. alignment information about both read pairs in adjoining rows). 
    # The alignment tool might sort them for you, but watch out for how the sorting was done. 
    # If they are sorted by coordinates (like with STAR), you will need to use samtools sort to re-sort them by 
    # read name before using as input in featureCounts. If you do not sort you BAM file by read name before using as input, 
    # featureCounts assumes that almost all the reads are not properly paired.

    # if [[ $read2 != "none" ]]
    # then
    #     samtools sort -n -O bam "$3" -o "${3}.tmp"
    #     mv "$3" "${3/bam/coord.bam}"
    #     mv "${3}.tmp" "$3"
    # fi
    # 
    # changed picard to queryname above, should work
    
    if [[ ! -z $cullfinks ]]
    then
        echo "Cufflinks started $4"
        
        #cufflinks requires coordinate sorted bam file
        samtools sort -@ $5 -o "${3/.bam/.coord.bam}" "$3"
        
        
        #Cufflinks
        if [[ $(basename $read2) == "none" ]]
        then
            cufflinks -q -p $5 -o "$1" -m $7 -g "$4" "${3/.bam/.coord.bam}"
        #-m is average fragment length - ie. for unpaired reads only
        else
            cufflinks -q -p $5 -o "$1" -g "$4" "${3/.bam/.coord.bam}"
        fi
        #CuffQuant to ref
        cuffquant -q -p $5 -o "$1" "$4" "$3"
        #echo "seqname	source	feature	start	end	score	strand	frame	attributes" > "${read_file}.transcripts.gtf"
        #grep exon transcripts.gtf >> "${read_file}.exon.transcripts.gtf"
        #rename files
        mv "${1}/abundances.cxb" "${3/.coord.bam/.abundances.cxb}"
        mv "${1}/genes.fpkm_tracking" "${3/.coord.bam/.genes.fpkm_tracking}"
        mv "${1}/isoforms.fpkm_tracking" "${3/.coord.bam/.isoforms.fpkm_tracking}"
        mv "${1}/skipped.gtf" "${3/.coord.bam/.skipped.gtf}"
        mv "${1}/transcripts.gtf" "${3/.coord.bam/.transcripts.gtf}"
        echo "Cufflinks completed"
        
        rm "${3/.bam/.coord.bam}"
    fi
    
    #get some stats such as number of mapped reads
    #this is outputted by star better but not by bowtie
    samtools flagstat "$3" > "${3/.bam/_flagstat.txt}"
    
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
    
    
    if [[ $strand == "reverse" ]]
    then
      stran_fc=2
      stran_qm="strand-specific-reverse"
    elif [[ $strand == "yes" ]]
    then
      stran_fc=1
      stran_qm="strand-specific-forward"
    else
      stran_fc=0
      stran_qm="non-strand-specific"
    fi
    
    if [[ $read2 != "none" ]]
    then
        fCount='-p' #this sets it to PE
    fi
    
    if [[ $6 == "B" ]]
    then
        htseq-count --type "gene" --idattr "Name" --order "name" --stranded="$strand" -a 5 --nonunique all -f bam "$3" "$4" > "${3/bam/HTSeq.counts}" #
        # or gene? - let user input type to count
        if [[ ! -z $feat ]]
        then
            featureCounts -F -d 30 -s "$stran_fc" -t "gene" -g "Name" -O -Q 5 --ignoreDup -T $5 -a "$4" -o "${3/.bam/.featCount.counts}" "$3" "$fCount"
        fi
        
    # elif [[ $8 == "miRNA" ]]
    # then
    # grep "miRNA" "$4" > "${4/.g*/.miRNA.gtf}"
    # htseq-count --order "pos" --stranded="$strand" -f bam "$3" "${4/.g*/.miRNA.gtf}" > "${3/.bam/.HTSeq.counts}"
    # featureCounts --ignoreDup -T $5 -a "$4" -o "${3/.bam/.featCount.counts}" "$3"
    else
        htseq-count --order "name" --stranded="$strand" -f bam "$3" "$4" > "${3/bam/HTSeq.counts}"
        if [[ ! -z $feat ]]
        then
            featureCounts -F -d 30 -s "$stran_fc" --ignoreDup -T $5 -a "$4" -o "${3/bam/featCount.counts}" "$3" "$fCount"
        fi
    fi
    echo "Counts completed"
    
    #can also do qualimap
    # export PATH=/users/bi/jlimberis/bin/qualimap_v2.2.1:$PATH
    if [[ ! -z $qualimap ]]
    then
        if [[ $(basename $read2) == "none" ]]
        then
            qualimap rnaseq -bam "$3" -gtf "$4" -outdir "${3/.bam/_qualimap}"
            qualimap comp-counts -bam "$3" -gtf "$4" -id "Name" -type "gene" -s -out "${3/.bam/_counts.html}"
        else
            qualimap rnaseq --paired --sorted -p "$stran_qm" -bam "$3" -gtf "$4" -outdir "${3/.bam/_qualimap}"
            qualimap comp-counts -bam "$3" -gtf "$4" -id "Name" -type "gene" -s -out "${3/.bam/_counts.html}" -p "$stran_qm" -pe
        fi
    fi
}


miRNAaln () {
    #miRNA alignment
    bowtie2 -p $1 --non-deterministic --very-sensitive -x "$2" -U ${3} | samtools view -@ $1 -Sb - > "$4"
    # --very-sensitive = -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    # or
    # bowtie2 -p "$threads" -D 20 -R 3 -N 1 -L 16 -i S,1,0.50 --non-deterministic -x "$g1" -U ${read1} | samtools view -@ $threads -Sb - > ${read1/.f*/.bam}
    # or
    # bowtie2 -p "$threads" -N 1 -L 16 --local -x "$g1" -U ${read1} | samtools view -@ $threads -Sb - > ${read1/.f*/.bam}
    #allow one mismatch for later SNP calling, seed length is 16
    
    #sort and index
    samtools sort -n -@ $1 "$4" -o "${4/.bam/.sorted.bam}"
    # samtools index "${4/.bam/.sorted.bam}"
    mv "${4/.bam/.sorted.bam}" "$4"
}

# not working
VaraintCall () {
    #GATK doesnt listen and eats ram so
    jav_ram=$(echo "scale=2; $ram*0.7" | bc)
    export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
    # if [ ! -f gatk ]; then
        # command -v gatk >/dev/null 2>&1 || { echo >&2 "I require gatk but it's not installed. Aborting."; exit 1; }
        # echo "$GATK not found! Canntor run SNP calling"
    # else
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
        gatk \
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
        gatk \
            -T RealignerTargetCreator \
            -R $1 \
            -I "${2}.split.bam" \
            -o "${2}.split.bam.list"
            #-known indels if available.vcf
            
        #Perform realignment of the target intervals
        gatk \
            -T IndelRealigner \
            -R $1 \
            -I "${2}.split.bam" \
            -targetIntervals "${2}.split.bam.list" \
            -o "${2}.tmp2.snps.bam"
            
        rm "${2}.split.bam"
          
          
          # Variant calling
        gatk \
            -T HaplotypeCaller \
            -R ${1} \
            -I "${2}.tmp2.snps.bam" \
            -dontUseSoftClippedBases \
            -o "${3}.vcf"
          
          # rm "${2}.tmp2.snps.bam"
        mv "${2}.tmp2.snps.bam" "${2/.bam/.snps.bam}"
          
          
          #Filter - we recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3
        gatk \
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
    # fi
    rm $(ls "${3}.split"*)
          
    #reset java mem
    jav_ram=$(echo "scale=2; $ram*0.8" | bc)
    export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
}






############
# pipeline #

if [[ ! -z $container ]]
then
    cd $container
    singularity pull "$container" library://semiquan7/default/rna_seq_pipeline
    exit 0
fi

if [[ ! -z $get_metrics ]]
then
    cd $get_metrics
    multiqc "$get_metrics" -n $(basename $get_metrics)
    exit 0
fi



# setup variables
declare_globals "$@"
# Script_dir_tmp=$(dirname "$0")
Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"
ram_def=$(expr $threads \* 2)
ram="${ram_in:-$ram_def}"
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-no}"
trim_min=16
# trim="${trim:-Y}" #Y|N
# keep_unpaired="${keep_unpaired:-N}" #Y|N
# vc="${vc:-n}"

# PATHs in singularity container
TRIM=/usr/bin/Trimmomatic-0.38/trimmomatic-0.38.jar
adapterSE=/usr/bin/Trimmomatic-0.38/adapters/universal.fa
adapterPE=/usr/bin/Trimmomatic-0.38/adapters/TruSeq2-PE.fa
PICARD=/usr/bin/picard.jar
# cut_adapt_seq="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TTACTATGCCGCTGGTGGCTCTAGATGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGA"
# GATK=/usr/bin/gatk-4.1.0.0/GenomeAnalysisTK.jar

out_dir="${out_dir:-read_dir}"
mkdir "${out_dir}"
name="${name:-${read1/.f*/}}"

read1="$read_dir/$read1"
read2="${read2:-none}"
read2="$read_dir/$read2"

mkdir "${out_dir}/${name}"
out_dir="${out_dir}/${name}"
back_dir=${PWD}
cd "$out_dir"

#check if programs installed
command -v cufflinks >/dev/null 2>&1 || { echo >&2 "I require cufflinks but it's not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo >&2 "I require bcftools but it's not installed. Aborting."; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v STAR >/dev/null 2>&1 || { echo >&2 "I require STAR but it's not installed. Aborting."; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo >&2 "I require htseq but it's not installed. Aborting."; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "I require python2 but it's not installed. Aborting."; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo >&2 "I require featureCounts but it's not installed. Aborting."; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require bowtie2 but it's not installed. Aborting."; exit 1; }
if [ ! -z $cut_adapt  ]; then command -v cutadapt >/dev/null 2>&1 || { echo >&2 "I require cutadapt but it's not installed. Aborting."; exit 1; }; fi

command -v python3 >/dev/null 2>&1 || { echo >&2 "I may require python3 but it's not installed."; }

if [ ! -f "$TRIM" ]; then echo "$TRIM not found!"; exit 1; fi
if [ ! -f "$PICARD" ]; then echo "$PICARD not found!"; exit 1; fi
if [ ! -z $vc  ]
then 
    command -v gatk >/dev/null 2>&1 || { echo >&2 "I require gatk but it's not installed. Aborting."; exit 1; }
fi

# get/set references
if [[ -z $g1 ]]
then 
    echo "No input genome supplied!"
    usage
    exit 1
else
    get_reference "$g1" "g1" "$gt1" "gt1"
fi

if [[ ! -z $g2 ]]
then
    get_reference "$g2" "g2" "$gt2" "gt2"
fi


# create index
if [ $t1 == "E" ] && [ -z $is_mi ]
then
    STAR_index "$threads" "$g1" "$gt1"
elif [ $t1 == "B" ] || [ ! -z $is_mi ] #if its miRNA or B then use bowtie
then
    BOWTIE_index "$g1" "$threads" "$gt1"
else
    echo "no type given for refernece 1, assuming eukaryotic"
    STAR_index "$threads" "$g1" "$gt1"
fi


if [[ ! -z $g2 ]]
then
    if [ $t2 == "E" ] && [ -z $is_mi ]
    then
        STAR_index "$threads" "$g2" "$gt2"
    elif [ $t2 == "B" ] || [ ! -z $is_mi ] #if its miRNA or B then use bowtie
    then
        BOWTIE_index "$g2" "$threads" "$gt2"
    else
        echo "no type given for refernece 2, assuming eukaryotic"
        STAR_index "$threads" "$g2" "$gt2"
    fi
fi


#fastqc and trim
if [[ $(basename $read2) == "none" ]]
then
    qc_trim_SE "$read1" "$out_dir" $ram $threads "$adapterSE" $trim_min
else
    qc_trim_PE "$read1" "$read2" "$out_dir" $ram $threads "$adapterPE" $trim_min
fi


#look for correct gtf for miRNA else make it
if [[ ! -z $is_mi ]]
then
    g_ext=".gtf"
    if [[ ! -e ${g1/.f*/.miRNA"$g_ext"} ]]
    then
        grep "miRNA" $gt1 > ${g1/.f*/.miRNA"$g_ext"}
        gt1=${g1/.f*/.miRNA"$g_ext"}
    fi
    if [[ ! -e ${g2/.f*/.miRNA"$g_ext"} ]]
    then
        if [[ $g2 != "none" ]]; then
              grep "miRNA" $gt2 > ${g2/.f*/.miRNA"$g_ext"}
              gt2=${g2/.f*/.miRNA"$g_ext"}
        fi
    fi
fi


#alignments
if [[ $(basename $read2) == "none" ]]
then
    #SE
    if [ $t1 == "B" ]
    then
        BOWTIE_alignerSE "$read1" $threads "$g1" "$out_dir" "$name" $ram
    elif [[ ! -z $is_mi ]]
    then
        miRNAaln $threads $g1 $read1 "${out_dir}/${name}.$(basename $g1).bam"
    else
        STAR_align $threads "$g1" "$read1" "$out_dir" "$name" $ram "$gt1"
    fi
else
    # PE
    if [[ $t1 == "B" ]]
    then
        BOWTIE_alignerPE "$read1" $threads "$g1" "$out_dir" "$name" $ram "$read2"
    elif [[ ! -z $is_mi ]] # miRNA will ony be SE
    then 
        echo "Cant process PE miRNA reads"
        exit 1
    else
        STAR_align $threads "$g1" "$read1" "$out_dir" "$name" $ram "$gt1" "$read2"
    fi
fi


bam_file="${out_dir}/${name}.$(printf $(basename $g1) | cut -f 1 -d '.').bam"
#this takes the first 10000 reads and calculates the read length
read_length=$(zcat $read1 | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')

if [[ -z "$only_care" ]]
then
    do_calcs $out_dir $g1 $bam_file $gt1 $threads $t1 $read_length
    if [[ ! -z $vc ]]; then
        VaraintCall "$g1" "$bam_file" "${out_dir}/${name}" "${name}"
    fi
fi


# unaligned
if [[ -e $g2 ]]
then
    # gen=$(basename $g2)
    # read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
    # read2_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate2.fastq.gz"
    #SE
    if [ $t2 == "B" ] && [ -z $is_mi ] #if its miRNA or B then use bowtie
    then
        BOWTIE_alignerSE "${read1_unaligned}" "$threads" "$g2" "$out_dir" "$name" "$ram"
    elif [[ ! -z  $is_mi ]]; then
        miRNAaln $threads $g2 ${read1_unaligned} "${out_dir}/${name}.$(basename $g2).bam"
    else
        STAR_align "$threads" "$g2" "${read1_unaligned}" "$out_dir" "$name" "$ram" "$gt2"
    fi
else
    # PE
    if [[ $t2 == "B" ]]
    then
        BOWTIE_alignerPE "$read1_unaligned" "$threads" "$g2" "$out_dir" "$name" "$ram" "$read2_unaligned"
    elif [[ ! -z  $is_mi ]] # miRNA will ony be SE
    then 
        echo "Cant process PE miRNA reads"
        exit 1
    else
        STAR_align "$threads" "$g2" "$read1_unaligned" "$out_dir" "$name" "$ram" "$gt2" "$read2_unaligned"
    fi
fi

if [[ -e $g2 ]]
then
    bam_file2="${out_dir}/${name}.$(printf $(basename $g2) | cut -f 1 -d '.').bam"
    #this takes the first 2500 reads and calculates the read length
    # read_length=$(zcat $read1_unaligned | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
    do_calcs "$out_dir" "$g2" "$bam_file2" "$gt2" $threads $t2 $read_length
    
    if [[ ! -z $vc ]]; then
        VaraintCall "$g2" "$bam_file2" "${out_dir}/${name}" "${name}"
    fi
fi

# Cleanup dirs





cd "$back_dir"





# cd "/home/lmbjas002/RNAseq_pipeline"; git pull
# declare -a samples=("C050_32183_CGTTGG" "C383_32190_TAAGAT" "C432_32197_GCCTGA" "C322_32184_TTCTGG" "C389_32191_AGGATC" "C440_32198_AGGTCG" "C353_32185_AACCGG" "C395_32192_ACTTAC" "C445_32199_CCTGCT" "C354_32186_CTGAGG" "C397_32193_CCGCAA" "C446_32200_TCATAC" "C378_32187_GTCTTC" "C402_32194_AATCCG" "C448_32201_GTCGCG" "C379_32188_ATTCTC" "C426_32195_GCAACG" "C380_32189_TGACTT" "C428_32196_TCAGAT")
# 
# R1="${samples[$SLURM_ARRAY_TASK_ID]}_read1.fastq.gz"
# R2="${samples[$SLURM_ARRAY_TASK_ID]}_read2.fastq.gz"
# name=${samples[$SLURM_ARRAY_TASK_ID]%%_*} #"${samples[0]:0:4}"
# readDir="/scratch/lmbjas002/CASS_sputum_crg"
# out_dir="/scratch/lmbjas002/CASS_sputumRun1"
# 
# singularity run /scratch/lmbjas002/rna_seq_pipeline_latest.sif bash /home/lmbjas002/RNAseq_pipeline/RNAseq_pipe_update.sh -t 8 \
#     --genome_reference1 "/scratch/lmbjas002/references/human/GCF_000001405.38_GRCh38.p12_genomic.fna" \
#     -g2 "/scratch/lmbjas002/references/tb/GCF_000195955.2_ASM19595v2_genomic.fna" \
#     --GTF_reference1 "/scratch/lmbjas002/references/human/GCF_000001405.38_GRCh38.p12_genomic.gtf" \
#     -gtf2 "/scratch/lmbjas002/references/tb/GCF_000195955.2_ASM19595v2_genomic.gtf" \
#     --Type_1 "E" \
#     -t2 "B" \
#     --read_dir "$readDir" \
#     --read1 "$R1" \
#     --read2 "$R2" \
#     -o "$out_dir" \
#     --name "$name" \
#     --feat_count \
#     --cufflinks \
#     --qualimap \
#     --strand "no" \
#     --trim \
#     --script_directory "/home/lmbjas002/RNAseq_pipeline" \
#     --only_care
# 
# 
# name2=${name}_only
# singularity run /scratch/lmbjas002/rna_seq_pipeline_latest.sif bash /home/lmbjas002/RNAseq_pipeline/RNAseq_pipe_update.sh -t 8 \
#     -g1 "/scratch/lmbjas002/references/tb/GCF_000195955.2_ASM19595v2_genomic.fna" \
#     -gtf1 "/scratch/lmbjas002/references/tb/GCF_000195955.2_ASM19595v2_genomic.gtf" \
#     --Type_1 "E" \
#     -t2 "B" \
#     --read_dir "$readDir" \
#     --read1 "$R1" \
#     --read2 "$R2" \
#     -o "$out_dir" \
#     --name "$name2" \
#     --feat_count \
#     --cufflinks \
#     --qualimap \
#     --strand "no" \
#     --trim \
#     --script_directory "/home/lmbjas002/RNAseq_pipeline"













#################
# random things #
# GFF to GTF 
# /Users/jdlim/bioinfomatics/cufflinks-2.2.1/gffread file.gff -T -o file.gtf
# perl ~/bin/gtf2bed.pl Homo_sapiens.GRCh38.86.gtf > Homo_sapiens.GRCh38.86.bed
# G1_bed=/users/bi/jlimberis/RNAseqData/ens_gen/Homo_sapiens.GRCh38.86.bed
# cd /users/bi/jlimberis/RNAseqData/ens_gen/trimmed/testing/T006
# read_distribution.py -i T006.miSub.fq.gz.sorted.bam -r $G1_bed
# 
# split reads into small and other RNA
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) <= 30 && length($2) >= 10) {print "@"$0} }' > ${i}.miRNA.fq
# bgzip ${i}.miRNA.fq
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) > 30) {print "@"$0} }' > ${i}.mRNA.fq
# bgzip ${i}.mRNA.fq
#################