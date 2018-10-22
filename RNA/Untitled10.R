#PBS -q UCTlong /
#PBS -l nodes=1:ppn=2:series600 /
#PBS -N CASS_rna_raw_count /
#PBS -V /
#PBS -M lmbjas002@myuct.ac.za /
#PBS -m ae /
cd "/researchdata/fhgfs/lmbjas002/CASS_RNAseq/Processed_out/bams"
for i in $(ls *.bam); do
/home/lmbjas002/bin/HTSeq/scripts/htseq-count -f bam $i /researchdata/fhgfs/lmbjas002/Server_Version/RNAseq/TBgenome/refs/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.28_renamed.gtf > ${i}.HTSeq.counts
done




cd "/researchdata/fhgfs/lmbjas002/CASS_RNAseq/processed/newly_processed"

for i in $(ls -d */); do
cd $i
name=${i/\//}
python2.6  /home/lmbjas002/bin/HTSeq/scripts/htseq-count -i ID -f bam "${name}_001.fastq.gz_accepted_hits.sorted.dedup.bam" "/researchdata/fhgfs/lmbjas002/Server_Version/RNAse
q/HumanGenome/Genome_ref/refs/GCF_000001405.34_GRCh38.p8_genomic.gff"> "${name}.HTSeq.counts"
cd ..
done



