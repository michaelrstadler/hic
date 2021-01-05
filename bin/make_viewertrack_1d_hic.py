"""

"""

from optparse import OptionParser
import sys
import re
import gzip
import io
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="paired alignment files, comma separated", metavar="FILE")
	
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")

	parser.add_option("-m", "--max_chr_size",
						dest="max_chr_size", default=1e8,
					  help="Max size to allocate for each chromosome. For many chromosomes (scaffolds), can get memory intensive if set high.")

	parser.add_option("-o", "--output_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds _binCounts_Xkb.txt")

	(options, args) = parser.parse_args()
	return options
                
options = parse_options()
bin_size = int(options.bin_size)
bin_counts = {}
max_bins = {}
filenames = options.filenames
files = filenames.split(',')
line_count = 0
max_chr_size = int(options.max_chr_size)
min_chr_size = 1e6
max_chr_bins = int(max_chr_size // bin_size)
min_chr_bins = int(min_chr_size // bin_size)

def add_read(chr_, bin_, bin_counts, max_chr_bins):
		if (chr_ not in bin_counts):
			bin_counts[chr_] = np.zeros(max_chr_bins)
		bin_counts[chr_][bin_] += 1

def update_max(chr_, bin_, max_bins):
	if (chr_ not in max_bins):
		max_bins[chr_] = bin_
	elif(bin_ > max_bins[chr_]):
		max_bins[chr_] = bin_

for f in files:
	print('Reading ' + f + '...')
	if (f[-2:] == 'gz'):
		file1 = gzip.open(f, 'rt')
	else:
		file1 = open(f, 'r')

	for line in file1:
		line_count = line_count + 1
		if (line_count % 10000000 == 0):
			print('.')
		line = line.rstrip()
		items = line.split()
		(chr1, Lmost1, chr2, Lmost2) = items[2], int(items[3]), items[5], int(items[6])
		bin1 = int(Lmost1 / bin_size)
		bin2 = int(Lmost2 / bin_size) 
		add_read(chr1,  bin1, bin_counts, max_chr_bins)
		add_read(chr2,  bin2, bin_counts, max_chr_bins)
		update_max(chr1, bin1, max_bins)
		update_max(chr2, bin2, max_bins)
	file1.close()
file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

outfilename = file_stem + '_1dBinCounts_' + str(bin_size) + '.txt'

outfile = open(outfilename, 'w')
for chr_ in bin_counts.keys():
	for bin_ in range(0, max_bins[chr_]):
		start = str(bin_ * bin_size)
		end = str((bin_ * bin_size) + bin_size - 1)
		val = str(bin_counts[chr_][bin_])
		outfile.write(chr_ + '\t' + start + '\t' + end + '\t' + val + '\n')
outfile.close()