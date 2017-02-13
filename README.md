Alignment of SE or PE Illumina sequenced RNA.
Script will align the reads to the supplied reference using:
STAR aligner for E(karyotyic) genomes 
Bowtie2 for P(prokaryotic) genomes

Each read file can first be aligned to one genome, then the unaligned sequences may be aligned to a second genome.
This is useful when comples, mixed clinical samples have be sequenced.

In the folder Deconvolute are scripts to remove sequences that may have arisen from the incorrect organism if the above is not used or if there are more than two expected orgnaisms in the sample.

Further unaligned read can be converted to fasta format and blasted.
#cat fastq_unalignes.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > fastq_unalignes..fa
#or using the nifty seqtk 

ToDo

SNPs can also be called from aligned reads using the Snp_calling_fromRNAseq.sh script. However, there number of aligned reads for many genes will mean this can only be used as a informal, exploratory guide.

differnetial expression and randomForest predictors can be made using the shiny app XXX

eQLT analysis ca be done using the shiny app XXX