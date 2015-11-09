__author__ = 'karlysindy'

import sys

alignfile = open(sys.argv[1], 'r')

OutFileName1 = "exactbar.bed"
OutFile1 = open(OutFileName1, 'w')

BarcodeArray=['GTTATGAAGG', 'ATCACTTAAG', 'ATAGCTCAGA', 'CGCCCTCGCA', 'ACCAAAAAAC', 'GACGGGGGTG', 'TCGAAACATA', 'GAGTCGTCTG', 'ATCTACCTGA', 'TCTGCTAGTT']

for line in alignfile:
    lines = line.strip('\n')
    linepart = line.split()
    bar = linepart[3]

    matching = [s for s in BarcodeArray if bar in s]

    if matching != []:

        OutFile1.write(line)
