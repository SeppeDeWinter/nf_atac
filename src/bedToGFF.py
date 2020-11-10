#!/usr/bin/env python

import sys
import argparse
import os.path

__author__ = 'kdavie'

""" Takes a bed file as input and converts it to a GFF format """

# Help and arguments parser
parser = argparse.ArgumentParser(description='Takes a bed file as input and\
 converts it to a GFF format')
parser.add_argument('inputBed', metavar='input', type=argparse.FileType('r'),
                    help='Input bed file', default=sys.stdin)
parser.add_argument('outputGFF', metavar='output', type=argparse.FileType('w'),
                    help='Output GFF path',
                    default=sys.stdout)
args = parser.parse_args()

# Check to see if bed is valid
inFileName, inFileExt = os.path.splitext(args.inputBed.name)
outFileName, outFileExt = os.path.splitext(args.outputGFF.name)

if inFileExt != '.bed':
    line = args.inputBed.readline().split()
    if not int(line[1]) > 0 or not int(line[2]) > 0:
        print('The input does not seem to be in a BED format. Either column\
              2 or 3 does not contain a number')
        exit()

# Sort the bed file by chromosome and then by start position
lines = [line.split() for line in args.inputBed]
lines.sort(key=lambda i: (i[0], int(i[1])))
regionNr = 1
newGFF = []
#BED is 0-based while GFF is 1-based, so have to add one to the coordinates
for line in lines:
    try:
        newLine = [line[0], 'Unk', 'Region', int(line[1]) + 1, int(line[2]) + 1, '.', '.', '.',
                   'region=', line[3]]
    except IndexError:
        newLine = [line[0], 'Unk', 'Region', int(line[1]) + 1, int(line[2]) + 1, '.', '.', '.',
                   'region=region_', regionNr]

    newGFF.append(newLine)
    regionNr += 1

for line in newGFF:
    args.outputGFF.write(str(line[0]) + '\t' + str(line[1]) + '\t' +
                         str(line[2]) + '\t' + str(line[3]) + '\t' +
                         str(line[4]) + '\t' + str(line[5]) + '\t' +
                         str(line[6]) + '\t' + str(line[7]) + '\t' +
                         str(line[8]) + str(line[9]) + '\n')
args.outputGFF.close()
