#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  -t|--threads
  -g1|--genome_reference1 = path to genome reference 1, if only a name is supplied the file will be downloaded from ncbi
  -g2|--genome_reference2 (optional) = path to genome reference 2, if only a name is supplied the file will be downloaded from ncbi
  
  Secondary alignment only works if first was E (will fix this sometime)
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
        is_mi="$2"
        ;;
        -c|--cufflinks) #Y or N
        cullfinks="$2"
        ;;
        -f|--feat_count) #Y or N
        feat="$2"
        ;;
        -q|--qualimap) #Y or N
        qualimap="$2"
        ;;
        -v|--variant_calling) #T or F
        vc="$2"
        ;;
        -s|--strand) #reverse
        strand="$2"
        ;;
        -tr|--trim) #Y|N
        trim="$2"
        ;;
        
    esac
        shift
    done
}


declare_globals "$@"
Script_dir=$(dirname "$0")
ram_def=$(expr $threads \* 2)
ram="${ram_in:-$ram_def}"
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-reverse}"
trim_min=10
trim="${trim:-Y}" #Y|N

echo $read2
