#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='gets kmer (eg. dimer) distribution for each position')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-k', required= True, help='k of kmer')
parser.add_argument('-n', required=  False, default=None, help='first n letters')
parser.add_argument('-l', required=  False, default=None, help='sequence length of interest')
parser.add_argument('--percentage', action='store_true', help = 'Write percentages instead if actual counts')

args = parser.parse_args()


if args.percentage:
    percentageFlag = True
else:
    percentageFlag = False

Fasta = fasta.fasta(args.i)
kmerAbundanceDict = Fasta.getKmerAbundance(int(args.k), args.n, args.l)


Fasta.writeKmerAbundanceTable(kmerAbundanceDict, args.o, percentageFlag)