from simplesam import Reader, Writer
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('r'), help="original SAM file")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="filtered SAM file")
parser.add_argument("-f", "--fasta", help="fasta file with coordinates in header")
args = parser.parse_args()
fastaFile = args.fasta
separator = '|'

def getReverseComplement(seq):
    def getComplementBase(base):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        if base in complement.keys():
            return complement[base]
        return base
    reverse_complement = "".join(getComplementBase(base) for base in reversed(seq))
    return reverse_complement

def decimalToBinary(n):
    return bin(n).replace("0b", "")

def getStrand(flag):
    binary = decimalToBinary(int(flag))
    if len(binary)>=5:
        if binary[-5] == "1":
            return "-"
    return "+"

fasta_in = open(fastaFile, 'r')
coordinates = set()
seqDict = {}
for line in fasta_in:
    if line.startswith('>'):
        header = line.strip()
        hl = header.split(separator)
        chromosome = hl[0]
        start = hl[1]
        end = hl[2]
        strand = hl[3]
        coordinateID = "_".join([start,end,strand])
        coordinates.add(coordinateID)
    else:
        seqDict[coordinateID] = line.strip()

in_file = args.input
in_sam = Reader(in_file)
out_file = args.output
out_sam = Writer(out_file, in_sam.header)

x = next(in_sam)
while True:
    if x.mapped:
        start = str(x.coords[0]-1)
        end = str(x.coords[-1])
        strand = getStrand(x.flag)
        item = "_".join([start,end,strand])
        if item in coordinates:
            if strand == '-':
                targetSequence = getReverseComplement(seqDict[item])
            else:
                targetSequence = seqDict[item]
            if x.seq == targetSequence:
                out_sam.write(x)
    try:
        x = next(in_sam)
    except StopIteration:
        print(x)
        break

out_sam.close()