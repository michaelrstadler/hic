'''
Takes a reduced Hi-C matrix file and implement "vanilla" normalization.

v2: outputs multiple bin sizes from a single input file. Hard-coded for 500 bp input data,
outputting data for 500, 1000, 2000, 4000 bp bins

Input needs to be for a single chromosome

No longer works with gzip

'''
from optparse import OptionParser
import sys
import re
import gzip
import os
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	(options, args) = parser.parse_args()
	return options

# function to retrieve the chromosome and maximum bin in the input file by peeking ahead to the end of the file.
def get_max_bin_chr(file):
	#This is somewhat complicated code that just finds the last line of the input file, returns the chromosome and bin number
	with open(file, 'rb') as f:
	    f.seek(-2, os.SEEK_END)
	    while f.read(1) != b'\n':
	        f.seek(-2, os.SEEK_CUR) 
	    last_line = (f.readline().decode())
	last_line_items = last_line.split()
	return((last_line_items[0], int(last_line_items[1])))            

# gets total counts for a bin of a new size (a combination of source bins)
def get_total_counts(bin_new, bins_combine, total_counts):
	total = 0
	for bin_orig in range(bin_new * bins_combine, (bin_new * bins_combine) + bins_combine):
		total = total + total_counts[bin_orig]
	return total

# gets total counts for a bin-bin linkage of a new size (a combination of source bins)
# essentially traces out a little square from the matrix
def get_bin_bin_counts(bin1_new, bin2_new, bins_combine, bin_bin_counts):
	count = 0
	for bin1_orig in range(bin1_new * bins_combine, (bin1_new * bins_combine) + bins_combine):
		for bin2_orig in range(bin2_new * bins_combine, (bin2_new * bins_combine) + bins_combine):
			count = count + bin_bin_counts[bin1_orig][bin2_orig]
	return count		

options = parse_options()
arb_constant = 10000000 #gets multiplied by the normalized coverage to avoid small number issues
width_limit = 1000 # won't print more bins than this in the "new" bin space. If fewer bins than 
# this are available for a given bin combination, will use the maximum number of bins possible.
bin_combinations = [1,2,4,8]
f = options.filename
(chromosome, max_bin) = get_max_bin_chr(f)
total_counts = np.zeros(max_bin + 1)
bin_bin_counts = np.zeros((max_bin+1, max_bin+1))
# handling of gzipped files
if (f[-2:] == 'gz'):
	file1 = gzip.open(f, 'rt')
else:
	file1 = open(f, 'r')

# prepare the file stem by removing extensions from input file
file_stem = re.sub('.txt', '', f)
file_stem = re.sub('.gz', '', file_stem)

#outfile = open(file_stem + '_VCnorm_bc' + str(bins_combine) + '.txt', 'w')
max_width = 0 #initialize max_width, will update when reading file

#read coverage file
for line in file1:
	line = line.rstrip()
	# the file is led by lines starting with '#' that contain the total counts recorded for each bin
	if (line[0] == '#'):
		line = line[1:]
		(chr, bin, count) = line.split('\t')
		total_counts[int(bin)] = int(float(count)) #int of '0.0' throws error...have to go through float for some reason
	# the lines that do not start with '#' contain the counts between bins
	else:
		(chr, bin1, bin2, count) = line.split('\t')
		bin1 = int(bin1)
		bin2 = int(bin2)
		bin_bin_counts[bin1][bin2] =  int(float(count)) #int of '0.0' throws error...have to go through float for some reason	
		width = abs(bin2 - bin1)
		if (width > max_width):
			max_width = width

#normalize and write file for various bin combinations

for bins_combine in bin_combinations:
	#create new file name, open this file.
	with open(re.sub('500', str(500 * bins_combine), file_stem) + '_VCnorm.txt', 'w') as outfile:
		# using consistent terminology where 'lnew' extension means bin numbers for this particular size (combination), whereas
		# the same variables without 'new' or with 'orig' are in the bin numbers of the input data
		max_bin_new = int(max_bin / bins_combine) - 1 #Including the last one can create problems if it isn't "full"--the program asks for bins that don't exist
		max_width_new = min(int(max_width / bins_combine), width_limit)
		for bin1 in range(0, max_bin_new + 1):
			for bin2 in range(max(0, (bin1 - max_width_new)), min(max_bin_new, bin1 + max_width_new + 1)):
				outfile.write(chromosome + '\t' + str(bin1)+ '\t' + str(bin2) + '\t')
				total_counts_bin1 = get_total_counts(bin1, bins_combine, total_counts)
				total_counts_bin2 = get_total_counts(bin2, bins_combine, total_counts)
				counts_bin1_bin2 = get_bin_bin_counts(bin1, bin2, bins_combine, bin_bin_counts)
				if (total_counts_bin1 != 0 and total_counts_bin2 != 0):
					count_norm = float(counts_bin1_bin2) / total_counts_bin1 / total_counts_bin2 * arb_constant
					outfile.write(str(round(count_norm, 2)) + '\n')
				else:
					outfile.write('NA\n')

file1.close()
outfile.close()