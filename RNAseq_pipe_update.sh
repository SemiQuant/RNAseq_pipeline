#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  -t|--threads
  -g1|--genome_reference1 = path to genome reference 1, if only a name is supplied the file will be downloaded from ncbi
  -g2|--genome_reference2 (optional) = path to genome reference 2, if only a name is supplied the file will be downloaded from ncbi
  "
}
    # echo -e 'Usage: ./RNA_processes.sh "1 - Input_paramaters.txt" "2 - threads" "3 - ram" \n
    # "4 - trim and clip adapt (Y|N)" "5 - is miRNA (Y|N)" "6 - stranded library (yes|no|reverse)" \n
    # "7 - cufflinks run in addition to HTseqCount (Y|N)" 8 - "SubRead FeatCount run in addition to HTseqCount (Y|N)" "9 - run qualimap (Y|N)" \n
    # "10 - SNP calling from seq data? (Y|N)" \n
    # "11 - genome1" "12 - GTF for genome 1" "13 - Type 1 (E=eukaryotic, B=bacterial)" \n
    # "14 - genome2" "15 - GTF for genome 2" "16 - Type 2 (E=eukaryotic, B=bacterial)" \n
    # if indexing a genome for the first time, this will require >30GB ram for a human genome\n'
    # 
    # echo -e "Input_paramaters.txt should be a comma seperated list conatining the following:\n
    # the directory where the read file(s) are, the output name, the output directory, the fastq file, \n
    # the second fastq file if PE reads - else leave blank (i.e file1,) \n"
    # exit 1
# }

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
    esac
        shift
    done
}

declare_globals "$@"
Script_dir=$(dirname "$0")


# supply path to genome files, you can leave blanks and only supply name if you wnat them downloaded, if they cannot be found they will be downloaded

get_reference () {
    mkdir "${Script_dir}/references" #wont overwrite so its ok
    if [[ ! -e $1 ]] #check if the file they supplied exists
    then
        local g1_f=$(basename $1)
        local g1_f=${g1_f%.f*}
        if [[ ! -e "${Script_dir}/references/${g1_f}.fasta" ]] #check if the file they supplied exists in the references folder
        then
            echo "Downloading reference genome ${g1_f}"
            curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${g1_f}&rettype=fasta" > "${Script_dir}/references/${g1_f}.fasta"
            export ${2}="${Script_dir}/references/${g1_f}.fasta"
        else
            export ${2}="${Script_dir}/${g1_f}.fasta"
            echo "Found reference genome file for ${g1_f}"
        fi
    else
        echo "Found reference genome file for $(basename $1)"
    fi
  
    if [[ ! -e $3 ]]
    then
        local gt1_f=$(basename $1)
        local gt1_f=${gt1_f%.f*}
        if [[ ! -e "${Script_dir}/references/${gt1_f}.gtf" ]]
        then
            curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${gt1_f}&rettype=gtf" > "${Script_dir}/references/${gt1_f}.gtf"
            export ${4}="${Script_dir}/references/${gt1_f}.gtf"
        else
            export ${4}="${Script_dir}/references/${gt1_f}.gtf"
            echo "Found annotations for genome file for ${gt1_f}"
        fi
    else
        echo "Found reference genome file for $(basename $1)"
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
        if [[ $3 == *.gz ]]
        then
            gunzip $3
        fi
  
    STAR \
      --runThreadN $1 \
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








############
# pipeline #
# set references
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
if [[ $t1 == "E" ]]
then
    STAR_index $threads $genome1 $G1
elif [[ $t1 == "B" ]]; then
    BOWTIE_index $genome1 $threads $G1
fi








# run
# /users/bi/jlimberis/CASS_RNAseq,C100,/users/bi/jlimberis/RNAseqData,C100_GTAGAG_HS374-375-376-merged_R1_001.fastq.gz,,
# /users/bi/jlimberis/testing/Homo_sapiens.GRCh38.dna.primary_assembly.fa,/users/bi/jlimberis/testing/GCF_000195955.2_ASM19595v2_genomic.fna,
# E,B,/users/bi/jlimberis/testing/Homo_sapiens.GRCh38.87.gtf,/users/bi/jlimberis/testing/GCF_000195955.2_ASM19595v2_genomic.gff

# singularity run ../RNAseq_pipe.sif bash ${PWD}/RNAseq_pipe_update.sh -t 8 \
# --genome_reference1 "/home/lmbjas002/RNAseq_pipeline/references/GCA_000001405.27_GRCh38.p12_genomic.fna" \
# -g2 "/home/lmbjas002/RNAseq_pipeline/references/GCF_000195955.2_ASM19595v2_genomic.fna" \
# --GTF_reference1 "/home/lmbjas002/RNAseq_pipeline/references/GCA_000001405.27_GRCh38.p12_genomic.gff" \
# -gtf2 "/home/lmbjas002/RNAseq_pipeline/references/GCF_000195955.2_ASM19595v2_genomic.gff" \
# --type1 "E" \
# -t2 "B"






