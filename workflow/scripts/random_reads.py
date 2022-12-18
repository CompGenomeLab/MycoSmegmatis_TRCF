import random
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="singleton fasta (.fa) or (.bed) file")
parser.add_argument("-b", "--bed", required=True, help="bed file to get the read lengths")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="bed output file")
args = parser.parse_args()
out = args.output

separator = '|'

bedIn = open(args.bed, 'r')
readLengthDict = {}
readLengthSet = set()
readCount = 0
for line in bedIn:
    ll = line.strip().split('\t')
    readLength = int(ll[2]) - int(ll[1])
    readLengthDict[readLength] = readLengthDict.get(readLength, 0) + 1
    readLengthSet.add(readLength)
    readCount += 1


bedIn.close()
masterDict = {}

if args.input.endswith('.fa'):
    fastaIn = open(args.input, 'r')
    for line in fastaIn:
        header = line.strip()[1:]
        seq = fastaIn.readline().strip()
        seqLen = len(seq)
        masterDict[seqLen] = masterDict.get(seqLen, [])
        masterDict[seqLen].append(header)
    fastaIn.close()
elif args.input.endswith('.bed'):
    for line in open(args.input, 'r'):
        hl = line.rstrip().split('\t')
        header = separator.join([hl[0], hl[1], hl[2], hl[5]])
        seqLen = int(hl[2]) - int(hl[1])
        masterDict[seqLen] = masterDict.get(seqLen, [])
        masterDict[seqLen].append(header)
else:
    assert "Unexpected file input"

# if args.input.endswith('.fa'):
#     for line in open(args.input, 'r'):
#         if line.startswith('>'):
#             singletonCount += 1
#             header = line.rstrip()
#         else:
#             seq = line.rstrip()
#             headerList.append(header[1:])
#             seqList.append(seq)
# elif args.input.endswith('.bed'):
#     for line in open(args.input, 'r'):
#         singletonCount += 1
#         hl = line.rstrip().split('\t')
#         header = separator.join([hl[0], hl[1], hl[2], hl[5]])
#         headerList.append(header)


selectedHeaders = []
for readLength in readLengthDict.keys():
    random.seed(int(readCount/readLength))
    headerList = masterDict[readLength]
    selectedHeaders += random.choices(headerList, k=readLengthDict[readLength])
        

for header in selectedHeaders:
    hl = header.split(separator)
    chromosome = hl[0]
    start = hl[1]
    end = hl[2]
    strand = hl[3]
    lineList = [chromosome, start, end, '.', '.', strand]
    out.write('\t'.join(lineList) + '\n')