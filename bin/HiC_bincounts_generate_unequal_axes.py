"""
 This script takes a paired alignment file and, assigns each end to two bins of two different
 supplied sizes. For example, a read at position 10,000,100 could be assigned to the 1 Mb bin 
 spanning 10M-11M and the 1 kb bin spanning 10.000 to 10.001 Mb. The result is a file in which
 each row represents a small bin and its linkages to the larger bins (columns), producing a
 matrix of unequal axes (say, 500 bp bins on the x axis and 50 kb bins on the y axis). The point
 is to assign each small bin a compartment based on its contact profile, avoiding the sparse
 matrix associated with using small bins on each axis.

 Prints a unique format. The file starts with a series of lines that start with a #.
 These are the total bin counts, used to do normalization, if desired, subsequently. After 
 that, the format is compressed:
 chromosome_bin    bin1  bin2  bin3 ...

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
	
	parser.add_option("-s", "--small_bin_size",
					   dest="small_bin_size", default=1000,
					  help="small bin size")
	
	parser.add_option("-l", "--large_bin_size",
					   dest="large_bin_size", default=25000,
					  help="large bin size")

	parser.add_option("-c", "--chromosome",
						dest="chromosome", default='none',
					  help="chromosome")
	
	parser.add_option("-o", "--outfile_stem",
						dest="outfile_stem", default='none',
					  help="output file stem. Adds ")

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
	if (bin in small_bin_totals):
		small_bin_totals[bin] = small_bin_totals[bin] + 1
	else:
		small_bin_totals[bin] = 1
	
	
options = parse_options()
small_bin_size = int(options.small_bin_size)
large_bin_size = int(options.large_bin_size)
selected_chromosome = options.chromosome
filenames = options.filenames
files = filenames.split(',')
bin_bin_counts =  {}
max_small_bin = 1
max_large_bin = 1
small_bin_totals = {}

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
			read1_small_bin = int(Lmost1 / small_bin_size)
			read1_large_bin = int(Lmost1 / large_bin_size)
			read2_small_bin = int(Lmost2 / small_bin_size)
			read2_large_bin = int(Lmost2 / large_bin_size)
			add_to_totals(read1_small_bin)
			if (read1_small_bin != read2_small_bin):
				add_to_totals(read2_small_bin) 

			if(read1_small_bin > max_small_bin): max_small_bin = read1_small_bin
			if(read2_small_bin > max_small_bin): max_small_bin = read2_small_bin
			if(read1_large_bin > max_large_bin): max_large_bin = read1_large_bin
			if(read2_large_bin > max_large_bin): max_large_bin = read2_large_bin

			Add_read(read1_small_bin, read2_large_bin)
			# Avoid double counting on diagonal. This will be triggered unless chromosome and bin are the same
			if (read1_small_bin != read2_small_bin):
					Add_read(read2_small_bin, read1_large_bin)

	infile.close()

file_stem = ''

if (options.outfile_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.outfile_stem

outfile = open(file_stem + '_unequalBinCounts_' + str(small_bin_size) + 'bp_' + str(large_bin_size) + 'bp_' + selected_chromosome + '.txt','w')

	 
print('done reading\n')

'''
for small_bin in range(0, max_small_bin + 1):
	outfile.write('#' + selected_chromosome + '\t' + str(small_bin) + '\t')
	if (small_bin in small_bin_totals):
		outfile.write(str(small_bin_totals[small_bin]) + '\n')
	else:
		outfile.write('0\n')
'''

for bin1 in range(0, max_small_bin + 1):
	outfile.write(selected_chromosome + '_' + str(bin1))
	if (bin1 in small_bin_totals):
		bin_total = small_bin_totals[bin1]
	else: bin_total = 1

	for bin2 in range(0, max_large_bin + 1):
		
		if (bin1 in bin_bin_counts):
			if(bin2 in bin_bin_counts[bin1]):
				bin_count_norm = float(bin_bin_counts[bin1][bin2]) / float(bin_total)
				outfile.write('\t' + str(bin_count_norm))
			else:
				outfile.write('\t0')
		else:
			outfile.write('\t0')
	outfile.write('\n')
outfile.close()