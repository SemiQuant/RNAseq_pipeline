#STAR aligner
threads=12
read_file=""

gen_dir="/researchdata/fhgfs/lmbjas002/Server_Version/RNAseq/HumanGenome/Genome_ref/refs"
cd $gen_dir
# Generating genome indexes
/home/lmbjas002/bin/STAR/bin/Linux_x86_64/STAR \
--runThreadN $threads \
--runMode genomeGenerate --genomeDir $gen_dir \
--genomeFastaFiles $gen_dir/GCF_000001405.34_GRCh38.p8_genomic.fa




#trim reads


#run alignment

/home/lmbjas002/bin/STAR/bin/Linux_x86_64/STAR \
--runThreadN $threads \
--genomeDir $gen_dir \
--sjdbGTFfile ${gen_dir}/GCF_000001405.34_GRCh38.p8_genomic.gtf \
--readFilesIn $read_file \
--readFilesCommand zcat
