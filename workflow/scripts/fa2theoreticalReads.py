import argparse
import sys

def getReverseComplement(seq):
    def getComplementBase(base):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        if base in complement.keys():
            return complement[base]
        return base
    reverse_complement = "".join(getComplementBase(base) for base in reversed(seq))
    return reverse_complement

def fa2theoreticalReads(filein, out, kmerSizeList, desiredCondition, desiredConditionReverse):
    separator = '|'
    seq = ''
    header = ''
    d = {}
    for line in filein:
        if line.strip().startswith('>'):
            if seq != '':
                d[header[1:]] = seq
                seq = ''
            header = line[1:].strip()
        else:
            seq += line.strip()
    d[header] = seq

    for header in d.keys():
        chromosome = header.split(' ')[0]
        seq = d[header]

        for k in kmerSizeList:
            for i in range(0,len(seq)-k):
                read = seq[i:i+k]
                if desiredCondition(read):
                    out.write('>' + chromosome + separator + str(i) + separator + str(i+k) + separator + '+' + '\n' + read.upper() + '\n')
                if desiredConditionReverse(read):
                    out.write('>' + chromosome + separator + str(i) + separator + str(i+k) + separator + '-' + '\n' + getReverseComplement(read).upper() + '\n')
    out.close()

def desiredCondition(seq):
    if seq[-4].upper() == 'T' and seq[-5].upper() == 'T':
        return True
    return False

def desiredConditionReverse(seq):
    if seq[3].upper() == 'A' and seq[4].upper() == 'A':
        return True
    return False

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('r'), help="Genome file in FASTA format")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="Theoretical reads in fasta format")
parser.add_argument('-k','--kmer', nargs='+', help='<Required> kmer length list', type=int, required=True)
args = parser.parse_args()

# kmerSizeList = [10, 11, 12, 13]
fa2theoreticalReads(args.input, args.output, args.kmer, desiredCondition, desiredConditionReverse)