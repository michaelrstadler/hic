'''
Converts WIG files (in 1, 2, and 4-column formats) to Hi-C track viewer
files. I believe it will handle different wig formats but I haven't 
pushed it much. It may break and need fixing with some files. 

Output format is chromosome [tab] bin [tab] value
'''
from optparse import OptionParser
import sys
import re
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="wig_file",
					  help="WIG file", metavar="FILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="out file", metavar="OUTFILE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="bin size", metavar="BINSIZE",
					  default=500)

	(options, args) = parser.parse_args()
	return options

def add_entry(chr_, bin_, value, bin_counts, max_bins, max_chr_bin):
	if (chr_ not in bin_counts):
		bin_counts[chr_] = np.zeros(max_chr_bin)
		max_bins[chr_] = 0
	bin_counts[chr_][bin_] += value
	if (bin_ > max_bins[chr_]):
		max_bins[chr_] = bin_

options = parse_options()
infile = open(options.wig_file,'r')
outfile = open(options.outfile, 'w')
bin_size = int(options.bin_size)
bin_counts = {}
max_bins = {}
max_chr_size = 1e8
max_chr_bin = int(max_chr_size // bin_size)

curr_chr = ''
curr_span = 0
curr_pos = 0
track_line_written = False

for line in infile:
	line = line.rstrip()
	if (line[0:5] == 'track'):
		pass
	elif(line[0:4] == 'vari'):
		(blank, chr, span) = line.split()
		curr_chr = chr[6:]
		if (not re.match('chr', curr_chr)):
			curr_chr = 'chr' + curr_chr
		curr_span = int(span[5:])
	elif(line[0:4] == 'fixe'):
		items = line.split()
		curr_chr = items[1]
		curr_chr = curr_chr[6:]
		curr_pos = 0
		curr_span = items[4]
		curr_span = int(curr_span[5:])
		if (not re.match('chr', curr_chr)):
			curr_chr = 'chr' + curr_chr
	elif(line[0] == '#'):
		pass
	else:
		line_values = line.split()
		if (len(line_values) == 4):
			(chr_, pos, pos2, value) = line.split()
		elif(len(line_values) == 2):
			(pos, value) = line.split()
			chr_ = curr_chr
		elif(len(line_values) == 1):
			value = line
			pos = curr_pos
			posR = posL + curr_span
			curr_pos = posR

		bin_ = int(pos) // bin_size
		value = float(value)
		add_entry(chr_, bin_, value, bin_counts, max_bins, max_chr_bin)

for chr_ in bin_counts.keys():
	for bin_ in range(0, max_bins[chr_]):
		val = bin_counts[chr_][bin_]
		if (val > 0):
			outfile.write(chr_ + '\t' + str(bin_) + '\t' + "{:.2f}".format(float(val)) + '\n')
infile.close()
outfile.close()




