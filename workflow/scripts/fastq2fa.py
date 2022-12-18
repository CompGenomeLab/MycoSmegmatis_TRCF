import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="FASTQ file")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="FASTQ file with removed duplicates")
args = parser.parse_args()
out = args.output


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
    out.write('>' + record['name'] + '\n' + record['seq'] + '\n')
