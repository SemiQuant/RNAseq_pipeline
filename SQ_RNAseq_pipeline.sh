#TODO
# get_reference () 
# write that this only takes gtfs and add in gff to gtf converter and notes
# add optinon to use cutadapt to remove set length from start or end
# also add option to supply adapter file ($cutAdapt declared but not used)
# add vaiables to pass other options to programmes
#     allow user to pass options to programes, like htseqCount
#   

# PATHs in singularity container
# Make the above an option ^ search for that line and you will know what I mean
#   
# also try add
# http://rseqc.sourceforge.net/
# https://bioconductor.org/packages/release/bioc/html/dupRadar.html
# http://smithlabresearch.org/software/preseq/
# http://hartleys.github.io/QoRTs/


# http://rseqc.sourceforge.net/#usage-information



# 
#   Notes
#     Think the REST for genomes has changed
#     lot of mixing of global and local variables, cleanup
#     add tmp dir
#     
#     rather, add option to remove it from the gff file
#     sed - '/rrl/d;/rrs/d;/rrf/d' test.txt ./infile
#     for tb its rrs, rrl, and rrf
#     
#     make functions delare local variable names that are descriptive
#     make it so you can choose to overwrite the star index
#     

usage () { #echo -e to read \n
  echo "
Usage Options
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Mode 1:
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -dl|--container = download the singularity container to this path and exit
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Mode 2:
@@@@@@@@@@@@@@@@';::;'@@@@@@@@@@@@@@@@@@    -mt|--get_metrics = supply a dir and get metrics for all analyses in that dir, all other options will be ignored if this is non-empyt
@@@@@@@@@@@@@:::::::::::::+@@@@@@@@@@@@@  Mode 3: 
@@@@@@@@@@@+:::::::::::::::::@@@@@@@@@@@    -ie|--index_exit = index genome 1 and exit, index type (Bowtie or Star) is based on Type_1 (Bacterial or Eukaryotic)
@@@@@@@@@@#::::::::::::::::::@@@@@@@@@@@  Mode 4: 
@@@@@@@@@@::::::'@@@@@#;::::@@@@@@@@@@@@    -t|--threads  
@@@@@@@@@@:::::'@@@@@@@@@@@+@@@@@@@@@@@@    -g1|--genome_reference1 = full path to genome reference (if only id then will downolad - currenlty not implemented due to REST  change)  
@@@@@@@@@@:::::'@@@@@@@@@@@@@@@@@@@@@@@@    -g2|--genome_reference2 (optional)   
@@@@@@@@@@::::::;@@@@@@@@@@@@@@@@@@@@@@@    -gtf1|--GTF_reference1  
@@@@@@@@@@@:::::::+@@@@@@@@@@@@@@@@@@@@@    -gtf2|--GTF_reference2  
@@@@@@@@@@@@::::::::;@@@@@@@@@@@@@@@@@@@    -t1|--Type_1 = E for eukaryotic or B for bacterial (defult is E)  
@@@@@@@@@@@@@+:::::::::@@@@@@@@@@@@@@@@@    -t2|--Type_2 = E for eukaryotic or B for bacterial (defult is E)  
@@@@@@@@@@@@@@@':::::::::'@@@@@@@@@@@@@@    -r|--ram  
@@@@@@@@@@@@@@+::::::::::::;@@@@@@@@@@@@    -rd|--read_dir  
@@@@@@@@@@@@:::::::::::::::::;@@@@@@@@@@    -r1|--read1  
@@@@@@@@@@:::::::::::''::::::::#@@@@@@@@    -r2|--read2 = the second fastq file if PE reads - else leave blank  
@@@@@@@@'::::::::+@@@@@@+:::::::;@@@@@@@    -o|--out_dir  
@@@@@@@::::::::@@@@@@@@@@@;:::::::@@@@@@    -n|--name  
@@@@@@:::::::@@@@@@@@@@@@@@@:::::::@@@@@    -m|--miRNA = Is this miRNA?  
@@@@@#::;:::@@@@@@@@@@@@@@@@@::::;:#@@@@    -c|--cufflinks = Run cufflinks?  
@@@@@::;:::@@@@@@@@@@@@@@@@@@@:::::;@@@@    -f|--feat_count = Run subread feature count?  
@@@@#::;::#@@@@@@@@@@@@@@@@@@@':::;;@@@@    -q|--qualimap = Run qualimap?  
@@@@:::;::@@@@@@@@@@@@@@@@@@@@@:::;;+@@@    -v|--variant_calling = (not implemented) Perform variant callineg? Currently not working
@@@@:::;::@@@@@@@@@@@@@@@@@@@@@::::;+@@@    -s|--strand = stranded library (yes|no|reverse)  
@@@@:::;::@@@@@@@@@@@@@@@@@@@@@:::;;+@@@    -tr|--trim = trim reads?
@@@@;::;::@@@@@@@@@@@@@@@@@@@@+:::;;@@@@    -sd|--script_directory  
@@@@@::;:::@@@@@@@@@@@@@@@@@@@:::::;@@@@    -fq|--fastQC = run fastqc?  
@@@@@::::::;@@@@@@@@@@@@@;;;;;:::::;@@@@    -sr|--shotRead = is read length short (like 50nt)?  
@@@@@@:::::::@@@@@@@@@@@;;;;;;;::::@@@@@    -sL|--SRlength = for making the index, if shotRead is on then this defult is 50  
@@@@@@@:::::::;@@@@@@@@@;;;;;;;:::@@@@@@    -rR|--remove_rRNA = remove rRNA from annotation file (not very robus, just deletes lines that say rRNA or ribosomal RNA), provide which GTF given it is (g1 or g2 or both)  
@@@@@@@@::::::::::'+##';,;;;;;;;:@@@@@@@    -rRm|--remove_rRNAmtb = remove rRNA from H37Rv annotation file, provide which GTF given it is (g1 or g2)  
@@@@@@@@@:::::::::::::::::;;;;;;;@@@@@@@    -k|--keep_unpaired
@@@@@@@@@@@:::::::::::::::::;;;;;;;@@@@@    -c2|--only_care = do you onlu care about genome 2?  
@@@@@@@@@@@@@@;::::::::::::+@@';;;;;+@@@    -mM|--multipleMet = picard multimet and rRNA met  
@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@#@@@@@    -s2|--star2 = basic 2 pass mode on star aligner?  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -ca|--cutAdapt = (not implemented) use cutadapt and pass options (such as '-u $nt_trim_from_start -m $minimum_length_keep -j $no_threads' )
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -hq|--htseq_qual = htseq-count quality cutoff (defult = 0)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -md2|--md5check1 = exit if read 1 md5 does not match this input 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -md2|--md5check2 = exit if read 2md5 does not match this input 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Mode 5:
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    -gtg|--GffToGtf = convert gff (path here) to gtf
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
"
}


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
        -mt|--get_metrics)
        get_metrics="$2"
        ;;
        -fq|--fastQC)
        fastQC="Y"
        ;;
        -sr|--shotRead)
        Sread="Y"
        ;;
        -sL|--SRlength)
        SRlen="50"
        ;;
        -rR|--remove_rRNA)
        rRNA="Y"
        ;;
        -rRm|--remove_rRNAMtb)
        rRNAmtb="$2"
        ;;
        -mM|--multipleMet)
        mMet="Y"
        ;;
        -s2|--star2)
        star2="Y"
        ;;
        -ie|--index_exit)
        ie="Y"
        ;;
        -ht|--htSeq)
        ht="Y"
        ;;
        -md1|--md5check1)
        md1="$2"
        ;;
        -md2|--md5check2)
        md2="$2"
        ;;
        -gtg|--GffToGtf)
        gtg="$2"
        ;;
        -ca|--cutAdapt)
        cutAdapt="$2"
        ;;
        -hq|--htseq_qual)
        htseq_qual="$2"
        ;;
    esac
        shift
    done
}


extract () {
    if [ -f $1 ] ; then
      case $1 in
        *.tar.bz2)   tar xjf $1     ;;
        *.tar.gz)    tar xzf $1     ;;
        *.bz2)       bunzip2 $1     ;;
        *.rar)       unrar e $1     ;;
        *.gz)        gunzip $1      ;;
        *.tar)       tar xf $1      ;;
        *.tbz2)      tar xjf $1     ;;
        *.tgz)       tar xzf $1     ;;
        *.zip)       unzip $1       ;;
        *.Z)         uncompress $1  ;;
        *.7z)        7z x $1        ;;
        *)  echo "'$1' is not zipped" 2>/dev/null ;;
         esac
    fi
}


check_md () {
    if md5sum -c <<<"$1  $2"
    then
        echo "md5s match $1"
    else
        echo "md5 Sums dont match!"
        echo "$1"
        exit 1
    fi
}


STAR_index () {
  # STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie" 2>&1 | tee -a "$log_file"
    if [[ ! -e "$(dirname $2)/chrLength.txt" ]] #check if indexed alread, not robust as require each to be in own folder..but so does star file naming..
    then
        if [ ! -z $5 ] || [ "$4" -lt 51 ]
        then
            if [ ! -z $6 ]
            then
                shot_read="--sjdbOverhang $(( $6 - 1 ))"
            else
                len="$4"
                shot_read="--sjdbOverhang $(( len - 1 ))"
            # --outFilterMatchNmin 0 --outFilterMismatchNmax 2
            fi
        fi
        
        # if [[ "${3##*.}" == "gff" ]]
        # then
        #     gff_used='--sjdbGTFtagExonParentTranscript Parent'
        # fi
        
        echo "Star indexing requires about 30GB ram for human genome, so if an error then check the log"
        STAR \
          --runThreadN "$1" \
          --runMode genomeGenerate \
          --genomeDir $(dirname $2) \
          --genomeFastaFiles "${2/.gz/}" \
          --sjdbGTFfile "${3/.gz/}" \
          --outFileNamePrefix ${2/.f*/} "$shot_read" #"$gff_used"
          # --sjdbOverhang read_length
          #this is the readlength of the RNA data - can get it from fastqs using fastqc or awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' fastq
          echo "STAR indexing complete on $2"
    else
        echo "Found STAR index for $2?"
    fi
    
    if [ ! -z $ie ]
    then
        exit
    fi
}


qc_trim_SE () {
    #FastQC pre
    if [[ ! -e "${3}/${1/.f*/_fastqc.zip}" ]]
    then
        if [[ ! -z $fastQC ]]
        then
            fastqc -t $5 "$1" -o "$3"
        fi
    fi
    
    #Trim Reads
    echo "trimming started $1"
    if [[ ! -z $trim ]]
    then
        read1="${3}/$(basename ${1/.f*/_forward.fq.gz})"
        if [[ -e "$read1" ]]
        then
            echo "Found ${1/.f*/}"
            export read1
        else
            java -jar "$TRIM" SE -phred33 \
              -threads $5 \
              "$1" \
              "${1/.f*/.trimmed.fq.gz}" \
              ILLUMINACLIP:"$6":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$7
      
            mv "${1/.f*/_forward_paired.fq.gz}" "$read1"
            
            #FastQC post
            if [[ ! -z $fastQC ]]
            then
                fastqc -t $5 "$read1" -o "$3"
            fi
        fi
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
      
            mv "${1/.f*/_forward_paired.fq.gz}" "$read1"
            mv "${2/.f*/_reverse_paired.fq.gz}" "$read2"
            
            if [[ ! -z $keep_unpaired ]]
            then
                #as we also want unpaired reads so..
                cat "${1/.f*/_forward_unpaired.fq.gz}" >> "$read1"
                cat "${2/.f*/_reverse_unpaired.fq.gz}" >> "$read2"
                
                # just to make sure as sometimes it lets one through if merging the paired and unpaired
                cutadapt --cores=$5 --quality-cutoff 10,10 --minimum-length $7 -o "${read1}.tmp.gz" "$read1"
                mv "${read1}.tmp.gz" "$read1"
                cutadapt --cores=$5  --quality-cutoff 10,10 --minimum-length $7 -o "${read2}.tmp.gz" "$read2"
                mv "${read2}.tmp.gz" "$read2"
            fi
            
            mv "${1/.f*/_forward_unpaired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" "${3}"
            
            #FastQC post
            if [[ ! -z $fastQC ]]
            then
                fastqc -t $5 "$read1" -o "$3"
                fastqc -t $5 "$read2" -o "$3"
            fi
        fi
        export read1
        export read2
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
  
  if [ ! -z $ie ]
  then
      exit
  fi
}


BOWTIE_alignerSE () {
    echo "BOWTIE alignment started $3"
    out_f="${4}/${5}.$(printf $(basename $3) | cut -f 1 -d '.').sam"
    
    if [[ -e "${out_f/.sam/.bam}" ]]
    then
        echo "Found ${out_f/.sam/.bam}"
        export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
    else
        bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -U "$1" -S "$out_f" --un-gz "${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
        
        export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz"
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
    local ref="${3/.f*/}"
    local gen=$(basename $ref)
    local out_f="${4}/${5}.$(printf $gen | cut -f 1 -d '.').sam"
    local out_f_bam="${out_f/sam/bam}"
    local R1="$1"
    local R2="$7"
    local thread=$2
    
    # local other_param=
    
    # -U <r> # Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the “standard in” or “stdin” filehandle.


    if [[ -e "$out_f_bam" ]]
    then
        echo "Found $out_f_bam"
        export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
        export read2_unaligned="${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"
    else
        # bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -1 "$1" -2 "$7" -S "$out_f" --un-gz ${4} --un-conc-gz ${4}
         bowtie2 \
            --dovetail \
            --local \
            --minins 0 \
            -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 \
            -p "$thread" -x "$ref" -1 "$R1" -2 "$R2" -S "$out_f" --un-gz "${4}" --un-conc-gz "${4}"

    #$(3 | cut -f 1 -d '.')
    # gen=$(basename $3)
    export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
    export read2_unaligned="${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"
       
    mv "un-conc-mate.1" "$read1_unaligned" 
    mv "un-conc-mate.2" "$read2_unaligned"
    
 
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
      O="$out_f_bam"\
      SORT_ORDER=queryname \
      VALIDATION_STRINGENCY=LENIENT

    rm "$out_f"
    
    fi
    echo "BOWTIE alignment completed"
}


STAR_align () {
    echo "Star alignment started"
    # STAR_align $threads "$g1" "$read1" "$out_dir" "$name" $ram "$gt1" "$read2" 2>&1 | tee -a "$log_file"
    local ref="${2/.f*/}"
    local ref_dir="$(dirname $2)"
    local gen=$(basename $ref)
    local thread=$1
    local R1="$3"
    local out_f="${4}/${5}.$(printf $(basename $2) | cut -f 1 -d '.').bam"
    local gtf="$7"
    local R2="$8"
    
    if [[ -e "$out_f" ]]
    then
        echo "Found ${out_f}"
    else
        if [[ "${gtf##*.}" == "gff" ]]
        then
            # if [[ ! -e "${gtf/gff/gtf}" ]]
            # then
                # gffread "$gtf" -T -o "${gtf/gff/gtf}"
            # fi
            # local gtf="${gtf/gff/gtf}"
             gff_used='--sjdbGTFtagExonParentTranscript Parent'
        fi
        
        
        if [[ ! -z $cufflinks ]]
        then
            CL='--outSAMstrandField intronMotif'
        fi
        
        if [ ! -z $star2 ]
        then
            TwoPass='--twopassMode Basic'
        fi

        #use two pass mode if intresited in novel jusctions..doubles runtime
        STAR \
          --runThreadN $thread \
          --genomeDir "$ref_dir" \
          --readFilesIn "$R1" "$R2" \
          --readFilesCommand zcat \
          --outFileNamePrefix "${4}/${5}" \
          --outSAMtype BAM Unsorted \
          --sjdbGTFfile "$gtf" \
          --outReadsUnmapped Fastx \
          --outSAMunmapped Within KeepPairs \
          --quantMode GeneCounts "$TwoPass" "$gff_used" "$CL" # "$shot_read" The counts coincide with those produced by htseq-count with default parameters. 
          # --outSAMunmapped
          
        # Unampped reads
        mv "${5}ReadsPerGene.out.tab" "${4}/${5}_ReadsPerGene.out.tab"
        #ovs this is only needed for PE but doesnt break anything
        mv "${4}/${5}Unmapped.out.mate1" "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
        bgzip "${4}/${5}_${gen}_Unmapped.out.mate1.fastq"
        mv "${4}/${5}Unmapped.out.mate2" "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"
        bgzip "${4}/${5}_${gen}_Unmapped.out.mate2.fastq"

        mv "${4}/${5}Aligned.out.bam" "$out_f"
        
    fi
    
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

    if [[ $strand == "reverse" ]]
    then
        stran_fc=2
        stran_qm="strand-specific-reverse"
        # LT=
    elif [[ $strand == "yes" ]]
    then
        stran_fc=1
        stran_qm="strand-specific-forward"
    else
        stran_fc=0
        stran_qm="non-strand-specific"
    fi
    
    if [[ ! -z $cullfinks ]]
    then
    
    # this is doing something to the bam file, like overwriting it
        echo "Cufflinks started $4"
        
        #cufflinks requires coordinate sorted bam file
        samtools sort -@ $5 -o "${3/bam/coord.bam}" "$3"
        
        #Cufflinks
        if [[ $(basename $read2) == "none" ]]
        then
            cufflinks -q -p $5 -o "$1" -m $7 -g "$4" "${3/bam/coord.bam}"
        #-m is average fragment length - ie. for unpaired reads only
        # –library-type "$LT" 
        else
            cufflinks -q -p $5 -o "$1" -g "$4" "${3/bam/coord.bam}"
        fi
        # CuffQuant to ref
        
        
        # has to be sam file??
        cuffquant --quiet --num-threads $5 --output-dir "$1" "$4" "${3/bam/coord.bam}"
        # echo "seqname	source	feature	start	end	score	strand	frame	attributes" > "${read_file}.transcripts.gtf"
        # grep exon transcripts.gtf >> "${read_file}.exon.transcripts.gtf"
        # rename files
        mv "${1}/abundances.cxb" "${3/coord.bam/abundances.cxb}"
        mv "${1}/genes.fpkm_tracking" "${3/coord.bam/genes.fpkm_tracking}"
        mv "${1}/isoforms.fpkm_tracking" "${3/coord.bam/isoforms.fpkm_tracking}"
        mv "${1}/skipped.gtf" "${3/coord.bam/skipped.gtf}"
        mv "${1}/transcripts.gtf" "${3/coord.bam/transcripts.gtf}"
        echo "Cufflinks completed"
        
        if [[ -z $mMet ]]
        then
            rm "${3/bam/coord.bam}"
        fi
    fi
    
    #get some stats such as number of mapped reads
    #this is outputted by star better but not by bowtie
    
    samtools flagstat "$3" > "${3/.bam/_flagstat.txt}"
    
    if [[ $6 == "E" ]]
    then
        local htseq_type; htseq_type="exon"
        local htseq_id; htseq_id="gene_id"
        local htseq_mode; htseq_mode="union"
        # local htseq_other; htseq_other="-a 5"
    else
        local htseq_type; htseq_type="gene"
        local htseq_id; htseq_id="Name"
        local htseq_mode; htseq_mode="union"
        local htseq_other; htseq_other="-a 5 --nonunique all"
        local fCount_bact; fCount_bact='-t "gene" -g "Name"'
    fi
    
    if [ ! -z $ht ]
    then
        echo "Started htseq-count $(basename $2)"
        # if [ ! -e "${4}.htseq.tmp.gff" ]
        # then
        #     sed '/exon-TRNF-1/d' "$4" > "${4}.htseq.tmp.gff"
        #     sed -i '/exon-RNR1-1/d' "{4}.htseq.tmp.gff"
        # fi
        htseq-count -a "$htseq_qual" --order "name" --type "$htseq_type" --idattr "$htseq_id" --mode "$htseq_mode" --stranded "$strand" $htseq_other -f bam "$3" "$4" > "${3/bam/HTSeq.counts}"
    fi
    
    
    if [[ $read2 != "none" ]]
    then
        local fCount
        fCount='-p' #this sets it to PE
        QMpaired="--paired"
    fi
    
    local gtf; gtf="$4"
    local fCount_other; fCount_other="-d 30 --ignoreDup"
    

    if [[ ! -z $feat ]]
    then
        echo "Started featureCounts $(basename $2)"
        featureCounts $fCount -F "GTF" -a "$gtf" -s "$stran_fc" -T $5 $fCount_other $fCount_bact -o "${3/bam/featCount.counts}" "$3"
        
        featureCounts $fCount -F "GTF" -a "$gtf" -s "$stran_fc" -g "gbkey" -T $5 -o "${3/.bam/_biotype.featureCounts.txt}" "$3"
        echo -e "$3\nBiotypes" > "${3/.bam/_biotype.featureCounts_out.txt}"
        cut -f 1,7 "${3/.bam/_biotype.featureCounts.txt}" | tail -n +3 >> "${3/.bam/_biotype.featureCounts_out.txt}"
    fi
    
    echo "Counts completed"
    
    #can also do qualimap
    # export PATH=/users/bi/jlimberis/bin/qualimap_v2.2.1:$PATH
    if [[ ! -z $qualimap ]]
    then
        #qualimap rnaseq only works with gtf file??
        # samtools sort -n -@ $5 -o "${3/bam/coord.bam}" "$3"
        # --sorted is giving issues.. have to let it do it? this take much more time
        echo "Started qualimap rnaseq $(basename $2)"
        qualimap rnaseq $QMpaired -p "$stran_qm" -bam "$3" -gtf "$gtf" -outdir "${3/.bam/_qualimap}" 1>/dev/null
          # --sorted
            
        echo "Started qualimap comp-counts $(basename $2)"
        if [ ! -e "${gtf}.qmap.tmp" ]
        then
            sed 's/exon/CDS/g' "$gtf" > "${gtf}.qmap.tmp.gtf"
        fi
        
        qualimap comp-counts -bam "$3" -gtf "${gtf}.qmap.tmp.gtf" -id "gene_id" -type "CDS" -s -out "${3/.bam/_Qualimap_counts.txt}" -p "$stran_qm" $QMpaired 1>/dev/null
            # --sorted
        # rm "${gtf}.qmap.tmp.gtf"
    fi
}


Multi_met_pic () {
    # $1 = gff
    # $2 = threads
    # $3 bam
    # $4 name
    # $5 ref
    # $6 strand
    
    #PICARD requires coordinate sorted bam file
    if [[ ! -e "${3/bam/coord.bam}" ]]
    then
        samtools sort -@ $2 -o "${3/bam/coord.bam}" "$3"
    fi
    
    java -jar "$PICARD" CollectMultipleMetrics \
        I="${3/bam/coord.bam}" \
        O="${4}_multiple_metrics" \
        R="$5"

    local gen=$(basename $1)
    local gen=${gen/.g*/}
    local gtf_1"=$1"
    
    
    if [[ ! -e "${gtf_1%.*}.rRNA.gtf" ]]
    then
        if [[ "${gtf_1##*.}" == "gff" ]]
        then
            cat "$gtf_1" | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > "${gtf_1%.*}.rRNA.gff" # All rRNAs 
            gffread "${gtf_1%.*}.rRNA.gff" -T -o "${gtf_1%.*}.rRNA.gtf"
        else
            cat "$gtf_1" | grep 'gbkey "rRNA"' | grep -v "ribosomal RNA protein" > "${gtf_1%.*}.rRNA.gtf"
        fi
        cat "${gtf_1%.*}.rRNA.gtf" | sed 's/\texon\t/\tgene\t/g' > tmp_rRNA # Add gene lines
        cat "${gtf_1%.*}.rRNA.gtf" | sed 's/\texon\t/\ttranscript\t/g' >> tmp_rRNA # Add transcript lines
        cat tmp_rRNA >> "${gtf_1%.*}.rRNA.gtf"
        rm tmp_rRNA
        sed -i -e 's/$/ gene_biotype \"rRNA\"; transcript_biotype \"rRNA\"; gene_source \"ncbi\"; transcript_source \"ncbi\";/' "${gtf_1%.*}.rRNA.gtf" # Add to end of each line
    fi
    
    if [[ ! -e "${gtf_1%.*}.refFlat.txt" ]]
    then
        gtfToGenePred -genePredExt "${gtf_1%.*}.rRNA.gtf" refFlat.txt.tmp
        paste <(cut -f 12 refFlat.txt.tmp) <(cut -f 1-10 refFlat.txt.tmp) > "${gtf_1%.*}.refFlat.txt" # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
        rm refFlat.txt.tmp
        
        if [[ "${gtf_1##*.}" == "gff" ]]
        then
            gffread "$gtf_1" -T -o "${gtf_1%.*}.all.rRNA.gtf"
        else
            cp "$gtf_1" "${gtf_1%.*}.all.rRNA.gtf" 
            # ln -s ?
        fi
        
       
        sed -i '/unknown_transcript_/d' "${gtf_1%.*}.all.rRNA.gtf"
        gtfToGenePred -genePredExt "${gtf_1%.*}.all.rRNA.gtf" refFlat2.txt.tmp
        paste <(cut -f 12 refFlat2.txt.tmp) <(cut -f 1-10 refFlat2.txt.tmp) > "${gtf_1%.*}.all.refFlat.txt" # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
        rm refFlat2.txt.tmp
    fi
    
    
    # Intervals for rRNA transcripts.
    samtools view -H "$3" > "${4}_rRNA.intervalListBody.txt"

    cat ${gtf_1%.*}.rRNA.gtf | awk '$3 == "transcript"' | \
      cut -f1,4,5,7,9 | \
      perl -lane '
          /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
          print join "\t", (@F[0,1,2,3], $1)
      ' | \
      sort -k1V -k2n -k3n \
      >> "${4}_rRNA.intervalListBody.txt"
      
    
    if [[ "$6" == "reverse" ]]
    then
        local strandP="SECOND_READ_TRANSCRIPTION_STRAND"
    elif [[ $strand == "yes" ]]
    then
        local strandP="FIRST_READ_TRANSCRIPTION_STRAND"
    else
        local strandP="NONE"
    fi
    
        
    java -jar "$PICARD" CollectRnaSeqMetrics \
        I="${3/bam/coord.bam}" \
        O="${4}_RNA_Metrics" \
        REF_FLAT="${gtf_1%.*}.all.refFlat.txt" \
        STRAND_SPECIFICITY="$strandP" \
        RIBOSOMAL_INTERVALS="${4}_rRNA.intervalListBody.txt" \
        CHART_OUTPUT="${4}.pdf"
        
        #ASSUME_SORTED=FALSE

    rm "${3/bam/coord.bam}"
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


# VaraintCall "$g1" "$bam_file" "${out_dir}/${name}" "${name}"

############
# pipeline #
# setup variables
declare_globals "$@"


if [[ ! -z $container ]]
then
    cd $container
    singularity pull library://semiquan7/default/rna_seq_pipeline
    exit 0
fi


if [[ ! -z $get_metrics ]]
then
    cd $get_metrics
    multiqc "$get_metrics" -n $(basename $get_metrics)
    exit 0
fi


if [[ ! -z $gtg ]]
then
    gffread "$gtg" -T -o "${gtg/gff/_SQ_RNApipelin.gtf}"
fi




# PATHs in singularity container
TRIM=/usr/bin/Trimmomatic-0.38/trimmomatic-0.38.jar
adapterSE=/usr/bin/Trimmomatic-0.38/adapters/universal.fa
adapterPE=/usr/bin/Trimmomatic-0.38/adapters/TruSeq2-PE.fa
PICARD=/usr/bin/picard.jar

# Script_dir_tmp=$(dirname "$0")
Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"

ram_def=$(expr $threads \* 2)
ram="${ram_in:-$ram_def}"
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-no}"
trim_min=16


#####################
# create index - copy paste of below but for index and exit
if [ ! -z $ie ]
then
    if [ $t1 == "E" ] && [ -z $is_mi ]
    then
        STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie"
    elif [ $t1 == "B" ] || [ ! -z $is_mi ] #if its miRNA or B then use bowtie
    then
        BOWTIE_index "$g1" "$threads" "$gt1"
    else
        echo "no type given for refernece 1, assuming eukaryotic"
        STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie"
    fi
fi
#####################




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

log_file="${out_dir}/${name}.RNAseq_pipeline_SemiQuant.log"
touch "$log_file"

echo "Command line input: $@" 2>&1 | tee -a "$log_file"


# set some defults

if [[ -z "$htseq_qual" ]]
then 
    htseq_qual=0
fi

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
if [ ! -z $cutAdapt  ]; then command -v cutadapt >/dev/null 2>&1 || { echo >&2 "I require cutadapt but it's not installed. Aborting."; exit 1; }; fi

command -v python3 >/dev/null 2>&1 || { echo >&2 "I may require python3 but it's not installed."; }

if [ ! -f "$TRIM" ]; then echo "$TRIM not found!"; exit 1; fi
if [ ! -f "$PICARD" ]; then echo "$PICARD not found!"; exit 1; fi


# create index
# get/set references
if [[ -z "$g1" ]]
then 
    echo "No input genome supplied!"
    usage
    exit 1
# else
#     get_reference "$g1" "g1" "$gt1" "gt1"
fi

#this takes the first 10000 reads and calculates the read length
read_length=$(zcat $read1 | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
echo "calculated read length is $read_length" 2>&1 | tee -a "$log_file"


# if input is gff then convert to gff
if [[ "${gt1##*.}" == "gff" ]]
then
    if [[ ! -e "${gtf1/gff/_SQ_RNApipelin.gtf}" ]]
    then
        gffread "$gtf1" -T -o "${gtf1/gff/_SQ_RNApipelin.gtf}"
        echo "Converted gff input to gtf, see file: ${gtf1/gff/_SQ_RNApipelin.gtf}"
    fi
    gtf1="${gtf1/gff/_SQ_RNApipelin.gtf}"
    export $gtf1
    echo "Using gtf file: ${gtf1/gff/_SQ_RNApipelin.gtf} instead of supplied $gtf1"
fi



# create index
if [ $t1 == "E" ] && [ -z $is_mi ]
then
    STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie" 2>&1 | tee -a "$log_file"
elif [ $t1 == "B" ] || [ ! -z $is_mi ] #if its miRNA or B then use bowtie
then
    BOWTIE_index "$g1" "$threads" "$gt1" 2>&1 | tee -a "$log_file" 
else
    echo "no type given for refernece 1, assuming eukaryotic"
    STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie" 2>&1 | tee -a "$log_file"
fi

if [[ ! -z $g2 ]]
then
    if [ $t2 == "E" ] && [ -z $is_mi ]
    then
        STAR_index "$threads" "$g2" "$gt2" "$read_length" "$Sread" "$SRlen" 2>&1 | tee -a "$log_file"
    elif [ $t2 == "B" ] || [ ! -z $is_mi ] #if its miRNA or B then use bowtie
    then
        BOWTIE_index "$g2" "$threads" "$gt2" 2>&1 | tee -a "$log_file"
    else
        echo "no type given for refernece 2, assuming eukaryotic"
        STAR_index "$threads" "$g1" "$gt1" "$read_length" "$Sread" "$SRlen" "$ie" 2>&1 | tee -a "$log_file"
    fi
fi



# check md sums
if [[ ! -z "$md1" ]]
then
    check_md "$md1" "$read1" 2>&1 | tee -a "$log_file"
fi

if [[ ! -z "$md2" ]]
then
    check_md "$md2" "$read2" 2>&1 | tee -a "$log_file"
fi



#fastqc and trim
if [[ $(basename $read2) == "none" ]]
then
    qc_trim_SE "$read1" "holder" "$out_dir" $ram $threads "$adapterPE" $trim_min 2>&1 | tee -a "$log_file"
else
    qc_trim_PE "$read1" "$read2" "$out_dir" $ram $threads "$adapterPE" $trim_min 2>&1 | tee -a "$log_file"
fi



#alignments
if [[ $(basename $read2) == "none" ]]
then
    #SE
    if [ $t1 == "B" ]
    then
        BOWTIE_alignerSE "$read1" $threads "$g1" "$out_dir" "$name" $ram 2>&1 | tee -a "$log_file"
    elif [[ ! -z $is_mi ]]
    then
        miRNAaln $threads $g1 $read1 "${out_dir}/${name}.$(basename $g1).bam" 2>&1 | tee -a "$log_file"
    else
        STAR_align $threads "$g1" "$read1" "$out_dir" "$name" $ram "$gt1" 2>&1 | tee -a "$log_file"
    fi
else
    # PE
    if [[ $t1 == "B" ]]
    then
        BOWTIE_alignerPE "$read1" $threads "$g1" "$out_dir" "$name" $ram "$read2" 2>&1 | tee -a "$log_file"
    elif [[ ! -z $is_mi ]] # miRNA will ony be SE
    then 
        echo "Cant process PE miRNA reads"
        exit 1
    else
        STAR_align $threads "$g1" "$read1" "$out_dir" "$name" $ram "$gt1" "$read2" 2>&1 | tee -a "$log_file"
    fi
fi


bam_file="${out_dir}/${name}.$(printf $(basename $g1) | cut -f 1 -d '.').bam"
if [[ -z "$only_care" ]]
then
    do_calcs $out_dir $g1 $bam_file $gt1 $threads "$t1" $read_length 2>&1 | tee -a "$log_file"

    if [[ ! -z $vc ]]; then
        VaraintCall "$g1" "$bam_file" "${out_dir}/${name}" "${name}" 2>&1 | tee -a "$log_file"
    fi
    
    if [[ ! -v $mMet ]]
    then
        Multi_met_pic "$gt1" $threads "$bam_file" "$name" "$g1" "$strand" 2>&1 | tee -a "$log_file"
    fi
fi


# unaligned
if [[ ! -e $g2 ]]
then
    cd "$back_dir"
    exit
fi


if [[ $read2 == "none" ]]
then
    # gen=$(basename $g2)
    # read1_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate1.fastq.gz"
    # read2_unaligned="${out_dir}/${name}_${gen}_Unmapped.out.mate2.fastq.gz"
    #SE
    if [ $t2 == "B" ] && [ -z $is_mi ] #if its miRNA or B then use bowtie
    then
        BOWTIE_alignerSE "$read1_unaligned" "$threads" "$g2" "$out_dir" "$name" "$ram" 2>&1 | tee -a "$log_file"
    elif [[ ! -z  $is_mi ]]; then
        miRNAaln $threads $g2 $read1_unaligned "${out_dir}/${name}.$(basename $g2).bam" 2>&1 | tee -a "$log_file"
    else
        STAR_align "$threads" "$g2" "${read1_unaligned}" "$out_dir" "$name" "$ram" "$gt2" 2>&1 | tee -a "$log_file"
    fi
else
    # PE
    if [[ $t2 == "B" ]]
    then
        BOWTIE_alignerPE "$read1_unaligned" "$threads" "$g2" "$out_dir" "$name" "$ram" "$read2_unaligned" 2>&1 | tee -a "$log_file"
    elif [[ ! -z  $is_mi ]] # miRNA will ony be SE
    then 
        echo "Cant process PE miRNA reads"
        exit 1
    else
        STAR_align "$threads" "$g2" "$read1_unaligned" "$out_dir" "$name" "$ram" "$gt2" "$read2_unaligned" 2>&1 | tee -a "$log_file"
    fi
fi




bam_file2="${out_dir}/${name}.$(printf $(basename $g2) | cut -f 1 -d '.').bam"
#this takes the first 2500 reads and calculates the read length
# read_length=$(zcat $read1_unaligned | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
do_calcs "$out_dir" "$g2" "$bam_file2" "$gt2" $threads $t2 $read_length 2>&1 | tee -a "$log_file"

if [[ ! -v $mMet ]]
then
    Multi_met_pic "$gt2" $threads "$bam_file2" "$name" "$g2" "$strand" 2>&1 | tee -a "$log_file"
fi


if [[ ! -z $vc ]]; then
    VaraintCall "$g2" "$bam_file2" "${out_dir}/${name}" "${name}" 2>&1 | tee -a "$log_file"
fi







#################
# random things #
# GFF to GTF 
# /Users/jdlim/bioinfomatics/cufflinks-2.2.1/gffread file.gff -T -o file.gtf
# perl ~/bin/gtf2bed.pl Homo_sapiens.GRCh38.86.gtf > Homo_sapiens.GRCh38.86.bed
# G1_bed=/users/bi/jlimberis/RNAseqData/ens_gen/Homo_sapiens.GRCh38.86.bed
# cd /users/bi/jlimberis/RNAseqData/ens_gen/trimmed/testing/T006
# read_distribution.py -i T006.miSub.fq.gz.sorted.bam -r $G1_bed

# split reads into small and other RNA
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) <= 30 && length($2) >= 10) {print "@"$0} }' > ${i}.miRNA.fq
# bgzip ${i}.miRNA.fq
# zcat  $i | awk 'BEGIN {RS="\n@";FS="\n"} {if (length($2) > 30) {print "@"$0} }' > ${i}.mRNA.fq
# bgzip ${i}.mRNA.fq
#################


#############
# GRAVEYARD #
# 
# #look for correct gtf for miRNA else make it
# if [[ ! -z $is_mi ]]
# then
#     g_ext=".gtf"
#     if [[ ! -e ${g1/.f*/.miRNA"$g_ext"} ]]
#     then
#         grep "miRNA" $gt1 > ${g1/.f*/.miRNA"$g_ext"}
#         export gt1=${g1/.f*/.miRNA"$g_ext"}
#     fi
#     if [[ ! -e ${g2/.f*/.miRNA"$g_ext"} ]]
#     then
#         if [[ $g2 != "none" ]]; then
#               grep "miRNA" $gt2 > ${g2/.f*/.miRNA"$g_ext"}
#               export gt2=${g2/.f*/.miRNA"$g_ext"}
#         fi
#     fi
# fi

# 
# 
# # remove rRNA
# if [[ -v $rRNA ]]
# then
#     if [[ $rRNA == "g1" ]]
#     then
#         sed '/rRNA/d;/ribosomal RNA/d;/ribosomal/d' "$gt1" > "$gt1_no_rRNA"
#         export gt1="$gt1_no_rRNA"
#     elif [[ $rRNA == "g2" ]]
#     then
#         sed '/rRNA/d;/ribosomal RNA/d;/ribosomal/d' "$gt2" > "$gt2_no_rRNA"
#         export gt2="$gt2_no_rRNA"
#     else
#         sed '/rRNA/d;/ribosomal RNA/d;/ribosomal/d' "$gt1" > "$gt1_no_rRNA"
#         export gt1="$gt1_no_rRNA"
#         sed '/rRNA/d;/ribosomal RNA/d;/ribosomal/d' "$gt2" > "$gt2_no_rRNA"
#         export gt2="$gt2_no_rRNA"
#     fi
# fi

# 
# if [[ $rRNAmtb == "g1" ]]
# then
#     sed '/rrl/d;/rrs/d;/rrf/d' "$gt1" > "$gt1_no_rRNA"
#     export gt1="$gt1_no_rRNA"
# fi
# 
# 
# if [[ $rRNAmtb == "g2" ]]
# then
#     sed '/rrl/d;/rrs/d;/rrf/d' "$gt2" > "$gt2_no_rRNA"
#     export gt2="$gt2_no_rRNA"
# fi
#
#
# not working
# VaraintCall () {
#     #GATK doesnt listen and eats ram so
#     jav_ram=$(echo "scale=2; $ram*0.7" | bc)
#     export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
#     # if [ ! -f gatk ]; then
#         # command -v gatk >/dev/null 2>&1 || { echo >&2 "I require gatk but it's not installed. Aborting."; exit 1; }
#         # echo "$GATK not found! Canntor run SNP calling"
#     # else
#         #Add read groups, sort, mark duplicates, and create index
#         java -jar $PICARD AddOrReplaceReadGroups \
#             I="$2" \
#             O="${2}.tmp.snps.bam" \
#             SO=coordinate \
#             RGID="id" RGLB="library" RGPL="ILLUMINA" RGPU="machine" RGSM="${4}"
#           
#           #check if dict exists
#         if [ ! -f "${1/.f*/.dict}" ]; then
#             java -jar $PICARD CreateSequenceDictionary \
#               R="$1"
#               O="${1/.f*/.dict}"
#         fi
#         #check if fai exists
#         if [ ! -f "${1}.fai" ]; then
#             samtools faidx "$1"
#         fi
#           
#         #index bamfile
#         samtools index "${2}.tmp.snps.bam"
#           
#           
#         # #do this if it compalins
#         # java -jar $PICARD ReorderSam \
#         #     I="${2}.tmp.snps.bam" \
#         #     O="${2}.tmp.snps.reordered.bam" \
#         #     R="$1" \
#         #     CREATE_INDEX=TRUE
#         # rm "${2}.tmp.snps.bam"
#         # mv "${2}.tmp.snps.reordered.bam" "${2}.tmp.snps.bam"
#         # mv "${2}.tmp.snps.reordered.bai" "${2}.tmp.snps.bam.bai"
#           
#           
#         # Split'N'Trim and reassign mapping qualities
#         gatk \
#             -T SplitNCigarReads \
#             -R $1 \
#             -I "${2}.tmp.snps.bam" \
#             -o "${2}.split.bam" \
#             -rf ReassignOneMappingQuality \
#             -RMQF 255 \
#             -RMQT 60 \
#             -U ALLOW_N_CIGAR_READS
#           
#         rm "${2}.tmp.snps.bam"
#           
#         java -jar $PICARD BuildBamIndex \
#             I="${2}.split.bam" \
#             VALIDATION_STRINGENCY= LENIENT
#             
#           
#         #Create a target list of intervals to be realigned with GATK
#         gatk \
#             -T RealignerTargetCreator \
#             -R $1 \
#             -I "${2}.split.bam" \
#             -o "${2}.split.bam.list"
#             #-known indels if available.vcf
#             
#         #Perform realignment of the target intervals
#         gatk \
#             -T IndelRealigner \
#             -R $1 \
#             -I "${2}.split.bam" \
#             -targetIntervals "${2}.split.bam.list" \
#             -o "${2}.tmp2.snps.bam"
#             
#         rm "${2}.split.bam"
#           
#           
#           # Variant calling
#         gatk \
#             -T HaplotypeCaller \
#             -R ${1} \
#             -I "${2}.tmp2.snps.bam" \
#             -dontUseSoftClippedBases \
#             -o "${3}.vcf"
#           
#           # rm "${2}.tmp2.snps.bam"
#         mv "${2}.tmp2.snps.bam" "${2/.bam/.snps.bam}"
#           
#           
#           #Filter - we recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3
#         gatk \
#             -T VariantFiltration \
#             -R "${1}" \
#             -V "${3}.vcf" \
#             -window 35 \
#             -cluster 3 \
#             -filterName "GATK_recomm" -filter "FS > 30.0 || QD < 2.0" \
#             -o "${3}_filtered.vcf"
#           
#         #get coverage
#         bedtools genomecov -ibam "${2/.bam/.snps.bam}" -bga > "${2/.bam/.bed}"
#         bgzip "${2/.bam/.bed}"
#     # fi
#     rm $(ls "${3}.split"*)
#           
#     #reset java mem
#     jav_ram=$(echo "scale=2; $ram*0.8" | bc)
#     export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
# }




# Multi_met_pic () {
#     # $1 = gff
#     # $2 = threads
#     # $3 bam
#     # $4 name
#     # $5 ref
#     # $6 strand
#     
#     #PICARD requires coordinate sorted bam file
#     if [[ ! -e "${3/bam/coord.bam}" ]]
#     then
#         samtools sort -@ $2 -o "${3/bam/coord.bam}" "$3"
#     fi
#     
#     java -jar "$PICARD" CollectMultipleMetrics \
#         I="${3/bam/coord.bam}" \
#         O="${4}_multiple_metrics" \
#         R="$5"
#       
#       
#       
#     local gen=$(basename $1)
#     local gen=${gen/.g*/}
#     
#     if [[ ! -e "${1/.g?f/_refFlat.txt}" ]]
#     then
#         if [[ "${1##*.}" == "gff" ]]
#         then
#             if [[ ! -e "${1/.gff/.gtf}" ]]
#             then
#                 gffread "$1" -T -o "${1/.gff/.gtf}"
#             fi
#         fi
#         gtfToGenePred -genePredExt "${1/.gff/.gtf}" "${gen}_refFlat.tmp.txt"
#         paste <(cut -f 12 "${gen}_refFlat.tmp.txt") <(cut -f 1-10 "${gen}_refFlat.tmp.txt") > "${1/.g?f/_refFlat.txt}"
#         rm "${gen}_refFlat.tmp.txt"
#     fi
#     
#     # for only rRNA
#     if [[ ! -e "${1/.g?f/_rRNA_refFlat.txt}" ]]
#     then
#         cat "$1" | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > "${1%.*}.rRNA.gff" # All rRNAs
#         gffread "${1%.*}.rRNA.gff" -T -o "${1%.*}.rRNA.gtf"
#         cat "${1%.*}.rRNA.gtf" | sed 's/\texon\t/\tgene\t/g' > "${1%.*}_tmp_rRNA" # Add gene lines
#         cat "${1%.*}.rRNA.gtf" | sed 's/\texon\t/\ttranscript\t/g' >> "${1%.*}_tmp_rRNA" # Add transcript lines
#         cat "${1%.*}_tmp_rRNA" >> "${1%.*}.rRNA.gtf"
#         rm "${1%.*}_tmp_rRNA"
#         sed -i -e 's/$/ gene_biotype \"rRNA\"; transcript_biotype \"rRNA\"; gene_source \"ncbi\"; transcript_source \"ncbi\";/' "${1%.*}.rRNA.gtf" # Add to end of each line
#         # Prepare refFlat from gtf for Picard
#         gtfToGenePred -genePredExt "${1%.*}.rRNA.gtf" "${1%.*}.rRNA.gtf.refFlat.txt.tmp"
#         paste <(cut -f 12 "${1%.*}.rRNA.gtf.refFlat.txt.tmp") <(cut -f 1-10 "${1%.*}.rRNA.gtf.refFlat.txt.tmp") > "${1/.g?f/_rRNA_refFlat.txt}" # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
#         rm "${1%.*}.rRNA.gtf.refFlat.txt.tmp"
#     fi
#     
#     # local ref_flat="${1/.g?f/_refFlat.txt}"
#     # local ribo_int=${4}_rRNA.intervalListBody.txt
#     local ref_flat="${1/.g?f/_rRNA_refFlat.txt}"
#     local ribo_int="${4}_r_rRNA.intervalListBody.txt"
#     
#     # Intervals for rRNA transcripts.
#     samtools view -H "${3/bam/coord.bam}" > "$ribo_int"
#     
#     cat "$ref_flat" | awk '$3 == "transcript"' | \
#       cut -f1,4,5,7,9 | \
#       perl -lane '
#           /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
#           print join "\t", (@F[0,1,2,3], $1)
#       ' | \
#       sort -k1V -k2n -k3n \
#       >> "$ribo_int"
#     
#     
#     if [[ "$6" == "reverse" ]]
#     then
#         local strandP="SECOND_READ_TRANSCRIPTION_STRAND"
#     elif [[ $strand == "yes" ]]
#     then
#         local strandP="FIRST_READ_TRANSCRIPTION_STRAND"
#     else
#         local strandP="NONE"
#     fi
#         
#     java -jar "$PICARD" CollectRnaSeqMetrics \
#         I="${3/bam/coord.bam}" \
#         O="${4}_RNA_Metrics" \
#         REF_FLAT="$ref_flat" \
#         STRAND_SPECIFICITY="$strandP" \
#         RIBOSOMAL_INTERVALS="$ribo_int" \
#         CHART_OUTPUT="${4}.pdf"
#         #ASSUME_SORTED=FALSE
# #"${1/.g?f/_rRNA_refFlat.txt}"
# # "${1/.g?f/_refFlat.txt}"
#     rm "${3/bam/coord.bam}"
# }
# 



