import sys
import argparse
import math

parser = argparse.ArgumentParser(description='Converts operon txt file of "Operon-mapper" to bed file.')
parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-c', '--chromosome', type=str, help='chromosome of a bacterium')
parser.add_argument('-o','--output', type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()
out = args.output
args.input.readline() # Skip the title line

operonList = []
i = -1
for line in args.input:
    ll = line.strip().split('\t')
    if len(ll)==1:
        operonName = ll[0]
        operonStart = math.inf
        operonEnd = -math.inf
        operonStrand = '.'

        firstGeneLL = args.input.readline().strip().split('\t')
        operonStart = min(int(firstGeneLL[3]), operonStart)
        operonEnd = max(int(firstGeneLL[4]), operonEnd)
        operonStrand = firstGeneLL[5]
        i+=1
        operonList.append({'name': operonName, 'start': operonStart, 'end': operonEnd, 'strand': operonStrand})
    else:
        operonStart = min(int(ll[3]), operonStart)
        operonEnd = max(int(ll[4]), operonEnd)
        operonList[i]['start'] = operonStart
        operonList[i]['end'] = operonEnd
        
for operon in operonList:
    lineList = [args.chromosome, str(operon['start']), str(operon['end']), operon['name'], '.', operon['strand']]
    out.write('\t'.join(lineList) + '\n')
