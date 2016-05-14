#!/bin/python
import argparse
from Bio import SeqIO
__author__ = 'Jason Limberis'
parser = argparse.ArgumentParser(description='This is a program that will create subsequence fasta files of a given length using a 1bp walker for the input genome')
parser.add_argument('-l','--lengths', help='Desired length range x:y',required = True)
parser.add_argument('-I', '--Genome_fasta', help='Input genome in fasta format (currently canot handle multi-fasta files',required = True)
parser.add_argument('-G','--genes', help='Bed file of gene ranges (0-indexed, bed file default)',required = True)
parser.add_argument('-O', '--Out_dir', help='Desired output directory',required = True)
parser.add_argument('-B', '--one_based', help='if your bed file is one based then set this flas "T"', required = False)
args = parser.parse_args()
lengths = args.lengths.split(':')
for seq_record in SeqIO.parse(args.Genome_fasta, "fasta"):
    name = seq_record.id
    ref = seq_record.seq
file_name = args.Out_dir + '/' + name + args.lengths + '.fasta'
start_len = int(lengths[0])
while start_len <= int(lengths[1]):
    length = start_len
    start_len += 1
    with open(args.genes) as genes:
        for line in genes:
            line = line.split('\t')
            name = line[0]
            start = int(line[1]) #- 2
            if args.one_based == "T":
                start = start -1
            gene_length = int(line[2])-int(line[1])
            for i in range(gene_length-length):
                end = start + length
                name_n = '>' + name + '.' + str(i) + '_' + str(length) + '\n'
                seq = str(ref[start:end]) + '\n'
#               file_name = Out_dir + name + '_' + str(length) + '.fasta'
                start = start + 1
                with open(file_name, 'a') as out_file:
                    out_file.write(name_n)
                    out_file.write(seq)
print 'Now run alignment_create_devonvolute_bed'
