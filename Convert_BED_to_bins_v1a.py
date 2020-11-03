"""
 
"""

from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="paired alignment FILE", metavar="FILE")
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")

	(options, args) = parser.parse_args()
	return options

options = parse_options()
bin_size = int(options.bin_size)
max_bin = {}
counts = {}

file1 = open(options.filename, 'r')

for line in file1:
	line = line.rstrip()
	if (line[0:3] == 'chr'):
		items = line.split('\t')
		chr = items[0][3:]
		pos = int(items[1])
		value = float(items[3])
		bin = int(int(pos) / bin_size)
		if(counts.has_key(chr)):
			if (counts[chr].has_key(bin)):
				counts[chr][bin] = counts[chr][bin] + value
			else:
				counts[chr][bin] = value
		else:
			counts[chr] = {}
			counts[chr][bin] = value
		
		if(max_bin.has_key(chr)):
			if (bin > max_bin[chr]):
				max_bin[chr] = bin
		else:
			max_bin[chr] = bin

file1.close()
file_stem = re.sub('.txt', '', options.filename)

chromosomes = ['X','2L','2R','3L','3R','4']

outfile = open(file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt','w')

for chr in chromosomes:
	for i in range(0,max_bin[chr]):
		if(counts[chr].has_key(i)):
			outfile.write('chr' + chr + '_' + str(i) + '\t' + str(counts[chr][i]) + '\n')
		else:
			outfile.write('0\n')
			
outfile.close()