"""
 This script takes a paired alignment file and, assigns each end to a bin (some chunk of 
 a chromosome defined by supplied bin size), and prints out the bin-bin counts for only
 contacts within some width (controlled by width variable, currently set to 500 bins) of
 the diagonal (distance between the bins). 

 Prints a unique format. The file starts with a series of lines that start with a #.
 These are the total bin counts, used to do normalization, if desired, subsequently. After 
 that, the format is compressed:
 chromosome    bin1    bin2    counts


 
 Some programming notes to check for: handling double counting of diagonal. 
 
"""

from optparse import OptionParser
import sys
import re
import gzip

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

	parser.add_option("-c", "--chromosome",
						dest="chromosome", default='none',
					  help="chromosome")
	
	parser.add_option("-s", "--file_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds diag_bin_counts bin size and chr")

	(options, args) = parser.parse_args()
	return options
            

def Add_read (bin1, bin2):
	if (bin1 in bin_bin_counts):
		if(bin2 in bin_bin_counts[bin1]):
			bin_bin_counts[bin1][bin2] = bin_bin_counts[bin1][bin2] + 1
		else:
			bin_bin_counts[bin1][bin2] = 1
	else:
		bin_bin_counts[bin1] = {}
		bin_bin_counts[bin1][bin2] = 1

def add_to_totals(bin):
	if (bin in bin_totals):
		bin_totals[bin] = bin_totals[bin] + 1
	else:
		bin_totals[bin] = 1

	
options = parse_options()
bin_size = int(options.bin_size)
selected_chromosome = options.chromosome
filenames = options.filenames
files = filenames.split(',')
bin_bin_counts =  {}
max_bin = 1
width = int(options.width)
bin_totals = {}

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
		if ((selected_chromosome == chr1 == chr2) or ('chr' + selected_chromosome == chr1 == chr2)):
			bin1 = int(Lmost1 / bin_size)
			bin2 = int(Lmost2 / bin_size)
			add_to_totals(bin1)
			if (bin1 != bin2):
				add_to_totals(bin2) 
			'''update_max(bin1)
			update_max(bin2)'''
			if (bin1 > max_bin):
				max_bin = bin1
			if (bin2 > max_bin):
				max_bin = bin2
			if (abs(bin1 - bin2) <= width):
				Add_read(bin1, bin2)
				# Avoid double counting on diagonal. This will be triggered unless chromosome and bin are the same
				if (bin1 != bin2):
					Add_read(bin2, bin1)

	infile.close()

file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

outfile = open(file_stem + '_CompressedBinCounts_' + str(bin_size) + 'bp_' + selected_chromosome + '.txt','w')

	 
print('done reading\n')

for bin in range(0, max_bin + 1):
	outfile.write('#' + selected_chromosome + '\t' + str(bin) + '\t')
	if (bin in bin_totals):
		outfile.write(str(bin_totals[bin]) + '\n')
	else:
		outfile.write('0\n')

for bin1 in range(0, max_bin + 1):
	for bin2 in range(max(0, bin1 - width), min(max_bin, bin1 + width)):
		outfile.write(selected_chromosome + '\t' + str(bin1) + '\t' + str(bin2) + '\t')
		if (bin1 in bin_bin_counts):
			if(bin2 in bin_bin_counts[bin1]):
				outfile.write(str(bin_bin_counts[bin1][bin2]))
			else:
				outfile.write('0')
		else:
			outfile.write('0')
		outfile.write('\n')
outfile.close()