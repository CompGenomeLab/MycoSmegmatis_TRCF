import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="sorted FASTQ file")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="lenght distribution file (tab-separated)")
parser.add_argument("-nozero", "--nozerooutput", help="lenght distribution file (tab-separated), zero excluded")
args = parser.parse_args()
out = args.output
nozero_out = open(args.nozerooutput, 'w')

def fastqGenerator(filename):
    with open(filename) as fh:
        while True:
            name = fh.readline().rstrip() # name line
            seq = fh.readline().rstrip() # read base sequence
            placeholder = fh.readline() # placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(name) == 0:
                break
            yield {'name': name, 'seq': seq, 'placeholder':placeholder, 'qual': qual}


lengthDict = {}
for record in fastqGenerator(args.input):
    l = len(record['seq'])
    lengthDict[l] = lengthDict.get(l,0) + 1

for l in sorted(lengthDict.keys()):
    out.write(f'{l}\t{lengthDict[l]}\n')
    if l != 0:
        nozero_out.write(f'{l}\t{lengthDict[l]}\n')
nozero_out.close()
