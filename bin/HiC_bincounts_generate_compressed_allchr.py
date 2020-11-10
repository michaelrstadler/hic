"""
 This script takes a paired alignment file and, assigns each end to a bin (some chunk of 
 a chromosome defined by supplied bin size), and prints out the bin-bin counts for only
 contacts within some width of the diagonal (distance between the bins). 

 Prints a unique format. The file starts with a series of lines that start with a #.
 These are the total bin counts, used to do normalization, if desired, subsequently. After 
 that, the format is compressed:
 chromosome    bin1    bin2    count
 
 Some programming notes to check for: handling double counting of diagonal. 
 
 Suggested use: use this script to generate 500 bp bin data, use downscreen normalization script to
 combine bins to generate larger bins if desired. 
"""

from optparse import OptionParser
import sys
import re
import gzip
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="paired alignment files, comma separated", metavar="FILE")
	
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
	
	parser.add_option("-w", "--width",
					   dest="width", default=1000,
					  help="width in bins from diagonal")

	
	parser.add_option("-s", "--file_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds diag_bin_counts bin size and chr")

	(options, args) = parser.parse_args()
	return options
            

def Add_read (chr, bin1, bin2, num_bins):
	if (chr not in bin_bin_counts):
		bin_bin_counts[chr] = np.zeros((num_bins, num_bins))
	bin_bin_counts[chr][bin1][bin2] = bin_bin_counts[chr][bin1][bin2] + 1
	

def add_to_totals(chr, bin, num_bins):
	if (chr in bin_totals):
		bin_totals[chr][bin] = bin_totals[chr][bin] + 1
	else:
		bin_totals[chr] = np.zeros(num_bins)

def update_max_bin(chr, bin1, bin2):
	if (chr in max_bin):
		return(max(bin1, bin2, max_bin[chr]))
	else:
		return(max(bin1, bin2))
	
options = parse_options()
bin_size = int(options.bin_size)
filenames = options.filenames
files = filenames.split(',')
width = int(options.width)
num_bins = int(1e5) #dummy variable for array creation...chromosome entrants must be <50 million bp
bin_bin_counts = {}
bin_totals = {}
max_bin = {}

line_count = 0

for f in files:
	if (f[-2:] == 'gz'):
		infile = gzip.open(f, 'rt')
	else:
		infile = open(f, 'r')

	for line in infile:
		line_count = line_count + 1
		if (line_count % 10000000 == 0):
			print('. ' + str(line_count / 10000000))
		line = line.rstrip()
		items = line.split()
		(chr1, Lmost1, chr2, Lmost2) = items[2], int(items[3]), items[5], int(items[6])
		if ((chr1 == chr2)):
			bin1 = int(Lmost1 / bin_size)
			bin2 = int(Lmost2 / bin_size)
			add_to_totals(chr1, bin1, num_bins)
			if (bin1 != bin2):
				add_to_totals(chr1, bin2, num_bins) 
			max_bin[chr1] = update_max_bin(chr1, bin1, bin2)
			if (abs(bin1 - bin2) <= width):
				Add_read(chr1, bin1, bin2, num_bins) # Avoid double counting on diagonal. This will be triggered unless chromosome and bin are the same
				if (bin1 != bin2):
					Add_read(chr1, bin2, bin1, num_bins)
	infile.close()

file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

print('done reading\n')

for chr in bin_totals.keys():
	with open(file_stem + '_CompressedBinCounts_' + str(bin_size) + 'bp_' + str(chr) + '.txt','w') as outfile:
		#print bin counts
		for bin in range(0, max_bin[chr] + 1):
			outfile.write('#' + chr + '\t' + str(bin) + '\t')
			outfile.write(str(bin_totals[chr][bin]) + '\n')

		#print bin_bin counts
		for bin1 in range(0, max_bin[chr] + 1):
			for bin2 in range(max(0, bin1 - width), min(max_bin[chr] + 1, bin1 + width + 1)):
				outfile.write(chr + '\t' + str(bin1) + '\t' + str(bin2) + '\t' + str(bin_bin_counts[chr][bin1][bin2]) + '\n')
