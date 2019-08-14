# Alignment of SE or PE Illumina sequenced RNA
Script will align the reads to the supplied reference using:
STAR aligner for E(karyotyic) genomes 
Bowtie2 for P(prokaryotic) genomes

Each read file can first be aligned to one genome, then the unaligned sequences may be aligned to a second genome.
This is useful when complex, mixed clinical samples have be sequenced.

Further unaligned read can be converted to fasta format and blasted.
```bash
cat fastq_unalignes.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > fastq_unalignes..fa
```
or using the nifty seqtk 

## Usage
### Get git
```shell
git clone https://github.com/SemiQuant/RNAseq_pipeline.git
```
***
### Get contianer
```shell
bash ./SQ_RNAseq_pipeline.sh --container "./"
```
***
### Index genome and exit
```shell
singularity run RNAseq_pipe.sif bash ./RNApipeline.sh \
  --index_exit \
  --threads 20 \
  --genome_reference1 "Path_to_genome.fna" \
  --GTF_reference1 "Path_to_genome.gtf" \
  --SRlength "sequencing_length" \
  --Type_1 "E"
```
***
### Run a single sample with two genomes (remove the last two three lines to run only one genome)
```shell
R1="Sample_R1_001.fastq.gz"
R2="Sample_R2_001.fastq.gz"
nme="SampleTest"

singularity run RNAseq_pipe.sif bash ./SQ_RNAseq_pipeline.sh \
  --threads 8 \
  --genome_reference1 "Path_to_genome.fna" \
  --GTF_reference1 "Path_to_genome.gtf" \
  --Type_1 "E" \
  --read_dir "Path_to_directory_containing_read_files" \
  --read1 "$R1" \
  --read2 "$R2" \
  --out_dir "Path_to_output_directory" \
  --name "$nme" \
  --feat_count \
  --qualimap \
  --strand "yes" \
  --trim \
  --script_directory "../RNAseq_pipeline" \
  --fastQC \
  --multipleMet \
  --star2 \
  --keep_unpaired \
  --Type_2 "B" \
  --genome_reference2 "Path_to_genome2.fna" \
  --GTF_reference2 "Path_to_genome2.gtf"
```
***
### Example wrapper slurm
```shell
#!/bin/sh
#SBATCH --account=USER
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=6:00:00
#SBATCH --job-name="RNAseqRun_1"
#SBATCH --mail-user=jason.limberis@uct.ac.za
#SBATCH --array=0-6

declare -a samples=("sample1" "sample2" "sample3" "sample4" "sample5" "sample6" "sample7")

R1="${samples[$SLURM_ARRAY_TASK_ID]}_R1_001.fastq.gz"
R2="${samples[$SLURM_ARRAY_TASK_ID]}_R2_001.fastq.gz"
nme=${samples[$SLURM_ARRAY_TASK_ID]%%_*}
#${samples[$SLURM_ARRAY_TASK_ID]%%_*} #"${samples[0]:0:4}"

Add the above here
```
***
### Example wrapper portable batch system
```shell
#PBS -q UCTlong /
#PBS -l nodes=1:ppn=8:series600 /
#PBS -N CASS_SEQ /
#PBS -M jason.limberis@uct.ac.za /
#PBS -m ae /
#PBS -t 0-6 /

declare -a samples=("sample1" "sample2" "sample3" "sample4" "sample5" "sample6" "sample7")

R1="${samples[${PBS_ARRAYID}]}_R1_001.fastq.gz"
R2="${samples[${PBS_ARRAYID}]}_R2_001.fastq.gz"
nme=${samples[${PBS_ARRAYID}]%%_*}
#"${samples[${PBS_ARRAYID}]}" #"${samples[0]:0:4}"

Add the above here
```
***
### Get metrics of all files (after they have all been processed)
```shell
singularity run RNAseq_pipe.sif bash ./SQ_RNAseq_pipeline.sh \
    --get_metrics "../runs"
```
***
### Convert gff to gtf
```shell
singularity run RNAseq_pipe.sif bash ./SQ_RNAseq_pipeline.sh \
    --GffToGtf "Path_to_GFF.gff"
```


### Note
If the genome has not been indexed, this will happen and requires ~30Gb of ram.
Index will be built using the length of sequencing reads from the supplied file.


---
## Available options

| Flag | Description |
|-------------|-------------|
  -dl --container | download the singularity container to this path and exit |  
  -t, --threads |
  -g1, --genome_reference1 | full path to genome reference |
  -gtf1, --GTF_reference1 |
  -t1, --Type_1 | E for eukaryotic or B for bacterial (defult is E) |
  -g2, --genome_reference2 | full path to genome reference |
  -gtf2, --GTF_reference2 |
  -t2, --Type_2 | E for eukaryotic or B for bacterial (defult is E) |
  -r, --ram |
  -rd, --read_dir |
  -r1, --read1 | Read 1 name |
  -r2, --read2 | Read 2 name (leave blank is not PE) |
  -o, --out_dir | Output directory |
  -n, --name | Name of sequence to append to output files |
  -c, --cufflinks | Run cufflinks? |
  -f, --feat_count | Run subread feature count? |
  -q, --qualimap | Run qualimap? |
  -s, --strand | stranded library (yes, no, reverse) |
  -tr, --trim | trim reads? |
  -sd, --script_directory |
  -fq, --fastQC | run fastqc? |
  -sr, --shotRead | is read length short (like 50nt)? |
  -sL, --SRlength | for making the index, if shotRead is on then this defult is 50 |
  -ht, --htSeq | do htseq counts, default is on for Bacterial genome |
  -mM, --multipleMet | get multiple metrics (picard) |
  -md1,--md5check1 | exit if md5 does not match this input for read 1 |
  -md2,--md5check2 | exit if md5 does not match this input for read 2 |
  -m, --miRNA | Is this miRNA? |
  -k, --keep_unpaired | Keep reads that dont align as fastq output? |
  -c2, --only_care | do you onlu care about genome 2? |
  -s2, --star2 | basic 2 pass mode on star aligner? |
  -ca, --cutAdapt | (not implemented) use cutadapt and pass options (such as '-u $nt_trim_from_start -m $minimum_length_keep -j $no_threads' ) |
  -hq, --htseq_qual | htseq-count quality cutoff (defult = 0) |
  -rRm, --remove_rRNAmtb | remove rRNA from H37Rv annotation file, provide which GTF given it is (g1 or g2) |
  -rR, --remove_rRNA | remove rRNA from annotation file (not very robus, just deletes lines that say rRNA or ribosomal RNA), provide which GTF given it is (g1 or g2 or both) |
  -v, --variant_calling | (not implemented) Perform variant callineg? Currently not working |
  *Single use options* |
  -mt, --get_metrics | supply a dir and get metrics for all analyses in that dir, all other options will be ignored if this is non-empyt |
  -ie, --index_exit | index genome and exit |
  -s2, --star2 | basic 2 pass mode on star aligner? |
  -ht, --htSeq | do htseq counts, alinger does this already 'on the fly' |
  -gtg, --GffToGtf | convert gff (path here) to gtf |
---

### These are not implemented int his script, contact Jason.Limberis@uct.ac.za for fancier things

| Flag|Description|
| -------------|-------------|
-t, --Type_1 | E for eukaryotic or B for bacterial |
-rR, --remove_rRNA | remove rRNA from annotation file  (not very robus, just deletes lines that say rRNA or ribosomal RNA), provide which GTF given it is (g1 or g2 or both)|
-rRm, --remove_rRNAmtb | remove rRNA from H37Rv annotation file|

## Important notes ong gtf/gff files
Some programs can only handle GTF or GFF files. Both are downloaded using the --downloadHuman and the scripts will look for the one not supplied otherwise create it if necessary.
[Here is a description of the formats](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/)
Also, in some instances, gtf are only used to store rRNA.

| Program|File|
|-------------|-------------|
  STAR | GTF/GFF |
  Cufflinks/cuffquant (dep)|GTF/GFF|
  htseq|GTF/GFF|
  featureCounts|GTF/GFF|
  qualimap|GTF|
  Picard|GTF/GFF?|


![](https://lh6.googleusercontent.com/proxy/zTRcW8Fc_7klEmoL3kXHkHeV323tDJ2wQGxS7f_LWzWDAOInVSUfnf0itVPvOPjg7ruPz2emKDFOBFS7jg=s0-d "lol")



# ToDo
---
1. Differnetial expression and randomForest predictors can be made using the shiny app [XXX](www.)

2. eQLT analysis can be done using the shiny app [XXX](www.)

3. Add [fastQscreen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) as an option


The code error checking isnt robust, I'd recommend doing getting the container, reference, then indexing the reference while doing an initial run. Then you can multi-run.

