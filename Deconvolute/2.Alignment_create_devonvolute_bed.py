import pysam
import argparse
import subprocess #for shell processes
import os
__author__ = "Jason Limberis"
parser = argparse.ArgumentParser(description = 'This program ')
parser.add_argument('-R', '--References', help = 'Bacterial or RNA references, seperated by ;', required = False, type = str)
parser.add_argument('-I', '--Input_file', help = 'Input multi fasta files from Make short reads', required = True, type = str)
parser.add_argument('-t', '--threads', help = 'number of threads to use, RAM will be 2*threads', required = True, type = int)
parser.add_argument('-r', '--Ref_euk', help = 'Eukaryotic references (introns), seperated by ;', required = False, type = str)
parser.add_argument('-M', '--MapQ', help = 'MapQ that will be used when aligning reads so as to filter out of bed file now. Default=10', required = False, type = int)

#>> still need to write the code for this
parser.add_argument('-S', '--Own_ref', help = 'The reference from which make the multifasta, this will bw used to remove any reads that align to more than one gene', required = False, type = str)
parser.add_argument('-G', '--Own_ref_genes', help = 'Own reference coding genes bed file', required = False)
#####
args = parser.parse_args()
ram = args.threads*2
if args.MapQ:
    mapQ = args.MapQ
else:
    mapQ = 10


#############
##Functions##
#############
def bed_from_sam(samIN, name):
    file = open(name, 'wa')
    pysam.view("-bS", "-o"+name, samIN)
    pysam.sort(name, name)
    Bamname = name+".bam"
    pysam.index(Bamname)
    #delete Sam file
    os.remove(samIN)
    bamfile = pysam.AlignmentFile(Bamname, "rb")
    for read in bamfile.fetch():
        if read.mapq > mapQ:
            #this may eliminate reads that aligned more than once, to cout how may use the XS flag
            line = "%s\t%i\t%i\n" % (str(read).split("\t")[0], read.reference_start, read.reference_end)
            file.write(line)
    file.close()
#################
##Functions end##
#################

#parser.add_argument('-O', '--Output_File', help = 'Output File Name', required = False, type = str)
#parser.add_argument('-B', '--Bed_File', help = 'To remove BED file', required = True)
nameIN = "ReadsToDeconvo.bed"
if args.References:
    bact_ref = args.References.split(';')
    for ref in bact_ref:
        #search for bowite2 references and make if not present
        try:
#>>       #make this detect the extension and use that, e.g .fa files
            with open(ref.replace(".fasta", ".1.bt2")) as index:
               pass
        except IOError as e:
            print "no bowtie2 index found, creating one"
            subprocess.call("bowtie2-build %s %s" %(ref, ref.replace(".fasta", "")), shell = True)
        ####
        ref = ref.replace(".fasta", "")
        print "processing %s" %ref
        out_file = "test"
        subprocess.call("bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p %s --very-sensitive -x %s -f %s -S %s" %(args.threads, ref, args.Input_file, out_file), shell = True)
        bed_from_sam(out_file, nameIN)

####check if this works
if args.Ref_euk:
    euk_ref = args.Ref_euk.split(';')
    for ref in euk_ref:
        #search for bowite2 references and make if not present
        try:
            with open(ref.replace(".fasta", ".1.bt2")) as index:
               pass
        except IOError as e:
            print "no bowtie2 index found, creating one"
            subprocess.call("bowtie2-build %s %s" %(ref, ref.replace(".fasta", "")), shell = True)
        ####
        #Align Reads to eukaryotic genome with introns and alternative splicing taken into account (Tophat2)
        subprocess.call("tophat2 -o %s -p %i --no-coverage-search %s %s" %(out_dir, threads, ref, args.Input_file), shell = True)
        out_file = ""
        bed_from_sam(out_file, nameIN)

#remove any duplicates that might have arisen due to alignment to more than one genome
#keep both files incase you want to see the areas of alignment confusion
lines_seen = set() # holds lines already seen
outfile = open(nameIN+"unique.bed", "w")
for line in open(nameIN, "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()

