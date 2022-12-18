import argparse
import sys


def fa2filteredFa(fastaFileName, out, length):
    filein = open(fastaFileName, 'r')
    while True:
        header = filein.readline().rstrip()
        seq = filein.readline().rstrip()
        if len(seq) == length:
            out.write(header + '\n' + seq + '\n')
        if len(header)==0:
            break


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="INPUT FASTA")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="OUTPUT FILTERED FASTA")
parser.add_argument("-l", "--length", type=int, help="LENGTH OF THE SEQUENCE")
args = parser.parse_args()

fa2filteredFa(args.input, args.output, args.length)