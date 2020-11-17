'''
Filters a paired mapping file to remove genomic (non Hi-C) reads. Removes 'in-out' (strand of Left end is +, Right end is -) reads that span
a distance less than the supplied distance cutoff.
'''
from optparse import OptionParser
import sys
import re
import random
import gzip
from random import randint
from random import shuffle
from numpy import percentile

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--infile", dest="infile",
					  help="input file: paired mapped", metavar="INFILE")
	parser.add_option("-d", "--distance", dest="distance",
					  help="distance cutoff", metavar="DISTANCE")
	(options, args) = parser.parse_args()
	return options



options = parse_options()
dist_cutoff = int(options.distance)

f1 = options.infile
if (f1[-2:] == 'gz'):
		infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

out_stem = re.sub('.gz', '', options.infile)
out_stem = re.sub('.txt', '', out_stem)
outfile = open(out_stem + '_genomicFiltered.txt', 'w')

for line in infile:
	line = line.rstrip()
	items = line.split()
	Lmost_strand = ''
	Rmost_strand = ''
	chr1 = items[2]
	chr2 = items[5]
	pos1 = int(items[3])
	pos2 = int(items[6])
	if (chr1 == chr2):
		size = abs(pos1 - pos2)
		if (size < dist_cutoff):
			if (pos1 < pos2):
				Lmost_strand = items[1]
				Rmost_strand = items[4]
			if (pos1 > pos2):
				Lmost_strand = items[4]
				Rmost_strand = items[1]
			
			if (Lmost_strand == '+' and Rmost_strand == '-'):
				continue
	outfile.write(line + '\n')

infile.close()
outfile.close()