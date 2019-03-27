sudo singularity build --sandbox RNAseq_pipe docker://ubuntu
sudo singularity shell --writable RNAseq_pipe
  apt-get update && apt-get install -y build-essential
  apt install zlib1g-dev
  apt install libssl-dev



# RNAseq_pipe
wget https://github.com/alexdobin/STAR/archive/2.7.0e.tar.gz
tar -xzf 2.7.0e.tar.gz
cd STAR-2.7.0e/source
make STAR



wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2


tar xjf samtools-1.9.tar.bz2
tar xjf bcftools-1.9.tar.bz2
tar xjf htslib-1.9.tar.bz2

apt-get install libncurses-dev
apt-get install libbz2-dev
apt-get install liblzma-dev
cd samtools-1.*
./configure
make
make install


apt-get install libcurl4-openssl-dev
cd bcftools-1.*
./configure
make
make install


apt-get install libssl-dev
cd htslib-1.*
./configure
make
make install




apt-get install openjdk-8-jdk  #or lighter use apt-get install openjdk-8-jre
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
apt-get install zip
unzip fastqc_v0.11.8.zip
cd /root/FastQC
chmod 755 fastqc
ln -s ${PWD}/fastqc /usr/local/bin/fastqc
fastqc_v0.11.8.zi


wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
rm Trimmomatic-0.38.zip


wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar

apt-get install python3.6
apt-get install python2.7
apt-get -y install python3-pip
pip3 install --user --upgrade cutadapt


apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
#apt install python-pip #pip install HTSeq

wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
tar -zxvf bedtools-2.28.0.tar.gz
cd bedtools2
make
make install


wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5/bowtie2-2.3.5-linux-x86_64.zip
unzip bowtie2-2.3.5-linux-x86_64.zip
rm bowtie2-2.3.5-linux-x86_64.zip
ln -s /root/bedtools2/bowtie2-2.3.5-linux-x86_64/bowtie2 /usr/local/bin/bowtie2

wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -xzf cufflinks-2.2.1.Linux_x86_64.tar.gz
rm cufflinks-2.2.1.Linux_x86_64.tar.gz
cp cufflinks-2.2.1.Linux_x86_64/* /usr/local/bin/

wget https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
unzip gatk-4.1.0.0.zip
rm gatk-4.1.0.0.zip


wget https://sourceforge.net/projects/subread/files/subread-1.6.4/subread-1.6.4-Linux-x86_64.tar.gz
tar -xzf subread-1.6.4-Linux-x86_64.tar.gz
rm subread-1.6.4-Linux-x86_64.tar.gz


pip install multiqc
apt-get install curl


sudo singularity build RNAseq_pipe.sif RNAseq_pipe






