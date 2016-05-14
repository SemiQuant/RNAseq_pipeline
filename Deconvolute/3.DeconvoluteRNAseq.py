#!/bin/python
#------------------------------------------------------------------------------------------------
#         FILE:  Deconvolute_RNA.py
# REQUIREMENTS:  Python 2.x and depencies: pysam; argparse
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jason Limberis
#      COMPANY:  UCT
#      VERSION:  0.2
#      CREATED:  03.10.15
#     REVISION:  ---
#------------------------------------------------------------------------------------------------
import pysam
import argparse
#import subprocess #for shell processes
__author__ = "Jason Limberis"
parser = argparse.ArgumentParser(description = 'This program takes in the sam file outputted by bowtie2 from MiSeq reads and deconvolutes the reads, leaving only high confidence TB RNA reads')
parser.add_argument('-I', '--Input-File', help = 'Input SAM file', required = True, type = str)
parser.add_argument('-O', '--Output-File', help = 'Output File Name', required = False, type = str)
parser.add_argument('-B', '--Bed-File', help = 'To remove BED file', required = True)
parser.add_argument('-b', '--Is_BAM_File', help = 'If input is a binary alignment file then set this to "y"', required = False)
#parser.add_argument('-G', '--Genes', help = 'Gene ranges BED file', required = False)

# Create a variable to hold the parsed arguments
args = parser.parse_args()

#setup variables
samIN = args.Input_File
bedFileIN = args.Bed_File
if len(args.Output_File) > 0:
    outname = args.Output_File
else:
    outname = samIN.replace(".sam", "_deconvoluted")
###################################
##Create dictionary from bed file##
###################################
#for each start position make a dictionary entry with a dictionary having the end position and returing the read name
#> sort the file to speed up indexing..?
remove = {}
with open(bedFileIN, 'r') as bed:
		for line in bed:
				bed_line = line.strip().split()
				start = int(bed_line[1]) #pysam is 0-based
				end = int(bed_line[2]) #pysam gives end position plus one for .reference_end
				if start in remove:
						#add to that dictionary
						remove[start].update({end : bed_line[0]})
				else:
						#make new entry
						remove.update({start : {end : bed_line[0]}})
#######
##End##
#######

#import sam file and convert it to and indexed bam file (nesessary for pysam to do what we need)
if args.Is_BAM_File == "y":
    Bamname = samIN
else:
    pysam.view("-bS", "-o"+outname, samIN)
    pysam.sort(outname, outname)
    Bamname = outname+".bam"
pysam.index(Bamname)

#import bam file
bamfile = pysam.AlignmentFile(Bamname, "rb")

#bam file to bed file
outFile_rem_name = outname+"_removed.bam"
outFile_keep_name = outname+"_kept.bam"
#set up tempelates for out files
outFile_rem = pysam.AlignmentFile(outFile_rem_name, "wb", template=bamfile)
outFile_keep = pysam.AlignmentFile(outFile_keep_name, "wb", template=bamfile)

reads_removed_count = 0
for read in bamfile.fetch(): #add until_eof=True in () to make it not remove unaligned reads
#    if (read.is_reverse) and  read.reference_end in remove and read.reference_start in remove[read.reference_end]:
#        outFile_rem.write(read)
#        reads_removed_count += 1
#    if read.mapq > 0:  ##not having until_eof=True already does this
      if read.reference_start in remove and read.reference_end in remove[read.reference_start]:
         outFile_rem.write(read)
         reads_removed_count += 1
           #need to count how many reads are removed for each gene and also how much this reduces the FRPKM
      else:
         outFile_keep.write(read)
print reads_removed_count
bamfile.close()
outFile_rem.close()
outFile_keep.close()

#need to check reads that aligned more than once and see if they aligned in different genes
#if they did then must remove them
#genes = arg.Genes
#for read in bamfile.fetch():
#    if read.mapq == 0:
#        print read.pos
