#!/bin/python
import argparse
__author__ = 'Jason Limberis'
parser = argparse.ArgumentParser(description='This is a program that will create subsequence fasta files of a given length using a 1bp walker for the H37Rv genome')
parser.add_argument('-l','--length', help='Desired length',required=True)
args = parser.parse_args()
length = int(args.length)
ref_f = open('/Users/jdlim/Desktop/RNA/TB/Python_trial/H37Rv.fasta.oneline.txt', 'r')
ref = ref_f.read()
with open('/Users/jdlim/Desktop/RNA/TB/Python_trial/RefSeq.cds_H37Rv.csv') as genes:
    for line in genes:
        line = line.split(";")
        name = line[0]
        start = int(line[1]) - 2
        gene_length = int(line[2])-int(line[1])
        for i in range(gene_length-length):
            start = start + 1
            end = start + length
            name_n = '>' + name + '.' + str(i) + '_' + str(length) + '\n'
            seq = ref[start:end] + '\n'
            file_name = '/Users/jdlim/Desktop/RNA/TB/Python_trial/Reads/' + name + '_' + str(length) + '.fasta'
            with open(file_name, 'a') as out_file:
                out_file.write(name_n)
                out_file.write(seq)
