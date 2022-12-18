import argparse
import sys

def fa2singletonList(fa):
    filein = open(fa, 'r')
    seq = ''
    seqList = []
    for line in filein:
        if line.strip().startswith('>'):
            if seq != '':
                seqList.append(seq)
                seq = ''
        else:
            seq += line.strip()
    seqList.append(seq)
    seqList.sort()

    previousSeq = ''
    multiCopySet = set()
    for seq in seqList:
        if seq == previousSeq:
            multiCopySet.add(seq)
        previousSeq = seq

    allSet = set(seqList)
    singletons = allSet.difference(multiCopySet)
    return singletons

    

def fa2getSingletons(fastaFileName, out, bedoutputfile):
    separator = '|'
    singletons = fa2singletonList(fastaFileName)
    filein = open(fastaFileName, 'r')
    bedout = open(bedoutputfile, 'w')
    readSet = set()
    positionDict = {}
    seenBefore = []
    i = 0
    for line in filein:
        if line.strip().startswith('>'):
            header = line.strip()[1:]
            i += 1
        else:
            sequence = line.strip()
            if sequence in singletons:
                out.write('>' + header + '\n' + sequence + '\n')
                hl = header.split(separator)
                chromosome = hl[0]
                start = hl[1]
                end = hl[2]
                strand = hl[3]
                bedout.write('\t'.join([chromosome, start, end, '.', '.', strand]) + '\n')
    bedout.close()


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Theoretical reads file in FASTA format")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="Singleton reads in fasta format")
parser.add_argument("-b", "--bedoutput", help="Output in bed format")
args = parser.parse_args()

fa2getSingletons(args.input, args.output, args.bedoutput)