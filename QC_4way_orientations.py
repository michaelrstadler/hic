'''
does a few things. First, it counts the total number of the 4 orientations of read pairs (left-most first, ++, +-, -+, --). 
This is Erez's in-in, in-out, out-in, out-out . Second, it counts the same for only reads that with distances less than 2000 bp, 
and prints out a file with the distances for the four types.
'''
from optparse import OptionParser
import sys
import re
import gzip
import random
from random import randint
from random import shuffle
from numpy import percentile

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--infile", dest="filename",
					  help="input file: paired mapped", metavar="INFILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile stem", metavar="OUTFILE")
	(options, args) = parser.parse_args()
	return options



options = parse_options()
pp_count = 0
pm_count = 0
mp_count = 0
mm_count = 0
plus_plus = []
plus_minus = []
minus_plus = []
minus_minus = []

f = options.filename
if (f[-2:] == 'gz'):
	infile = gzip.open(f, 'rt')
else:
	infile = open(options.filename,'r')

for line in infile:
	items = line.split()
	Lmost_strand = ''
	Rmost_strand = ''
	chr1 = items[2]
	chr2 = items[5]
	pos1 = int(items[3])
	pos2 = int(items[6])
	if (chr1 == chr2):
		size = abs(pos1 - pos2)
		if (pos1 < pos2):
			Lmost_strand = items[1]
			Rmost_strand = items[4]
		if (pos1 > pos2):
			Lmost_strand = items[4]
			Rmost_strand = items[1]
		
		if (Lmost_strand == '+' and Rmost_strand == '+'):
			if(size < 2000): plus_plus.append(size)
			pp_count += 1
		if (Lmost_strand == '+' and Rmost_strand == '-'):
			if(size < 2000): plus_minus.append(size)
			pm_count += 1
		if (Lmost_strand == '-' and Rmost_strand == '+'):
			if(size < 2000): minus_plus.append(size)
			mp_count += 1
		if (Lmost_strand == '-' and Rmost_strand == '-'):
			if(size < 2000): minus_minus.append(size)
			mm_count += 1

infile.close()

max_pp = len(plus_plus)
max_pm = len(plus_minus)
max_mp = len(minus_plus)
max_mm = len(minus_minus)

max_entry = max(max_pp, max_pm, max_mp, max_mm)


print('2000 counts: ' + str(max_pp) + '\t' + str(max_pm) + '\t' + str(max_mp) + '\t' + str(max_mm) + '\n')
print('Total counts: ' + str(pp_count) + '\t' + str(pm_count) + '\t' + str(mp_count) + '\t' + str(mm_count) + '\n')
