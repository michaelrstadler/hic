"""
 This script takes a paired alignment file and, assigns each end to a bin (some chunk of 
 a chromosome defined by supplied bin size), and then prints out a matrix of bins vs. bins
 for all chromosomes.
 
 Some programming notes to check for: handling double counting of diagonal. 
 
 v2: does bin counting for each chromosome at once, automatically assigns file names. 
 v2.3: Prints bin names in the form of chrX_3 (3rd bin on X).
 v2.4: Allows selection of a single chromosome for higher-res work, multiple input files
"""

from optparse import OptionParser
import sys
import re
import gzip
import io

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="paired alignment files, comma separated", metavar="FILE")
	
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
					  
	parser.add_option("-c", "--chromosome",
						dest="chromosome", default='none',
					  help="chromosome")
	
	parser.add_option("-s", "--file_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds _binCounts_Xkb.txt")

	(options, args) = parser.parse_args()
	return options

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
            

def Add_read(chr1, bin1, chr2, bin2):
	if (chr1 in bin_bin_counts):
		if (bin1 in bin_bin_counts[chr1]):
			if(chr2 in bin_bin_counts[chr1][bin1]):
				if(bin2 in bin_bin_counts[chr1][bin1][chr2]):
					bin_bin_counts[chr1][bin1][chr2][bin2] = bin_bin_counts[chr1][bin1][chr2][bin2] + 1
					return
	#if the entry exists, it will have returned. If any of the ifs fails, we add it:
	bin_bin_counts[chr1][bin1][chr2][bin2] = 1

def update_max(chr, bin):
	if (chr in max_bin):
		if (bin > max_bin[chr]):
			max_bin[chr] = bin
	else:
		max_bin[chr] = bin

options = parse_options()
bin_size = int(options.bin_size)
selected_chromosome = options.chromosome
filenames = options.filenames
files = filenames.split(',')
bin_bin_counts =  AutoVivification()
max_bin = {}
line_count = 0

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
		if (selected_chromosome != 'none'):
			if (chr1 != selected_chromosome or chr2 != selected_chromosome):
				continue
		bin1 = int(Lmost1 / bin_size)
		bin2 = int(Lmost2 / bin_size) 
		update_max(chr1, bin1)
		update_max(chr2, bin2)
		Add_read(chr1, bin1, chr2, bin2)
		# Avoid double counting on diagonal. This will be triggered unless chromosome and bin are the same
		if (chr1 != chr1 or bin1 != bin2):
			Add_read(chr2, bin2, chr1, bin1)

	file1.close()

file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

chromosomes = []
outfile = ''
if (selected_chromosome == 'none'):
	chromosomes = ['X','2L','2R','3L','3R','4']
	outfile = open(file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt','w')
else:
	chromosomes = [selected_chromosome]
	outfile = open(file_stem + '_binCounts_' + str(bin_size) + 'bp_chr' + selected_chromosome + '.txt','w')



#print column labels
for chr3 in chromosomes:
	for bin3 in range(0, max_bin[chr3] + 1):
		outfile.write('\tchr' + chr3 + '_' + str(bin3))
outfile.write('\n')
	 

for chr1 in chromosomes:
	for bin1 in range(0, max_bin[chr1] + 1):
		outfile.write('chr' + chr1 + '_' + str(bin1) + '\t')
		for chr2 in chromosomes:	
			for bin2 in range(0, max_bin[chr2] + 1):
				count = bin_bin_counts[chr1][bin1][chr2][bin2]
				# Slightly janky, but with Autovivification, if it's empty it makes a new empty list and returns this.
				# Can do simple int test to distinguish
				if isinstance(count,int):
					outfile.write(str(count))
				else:
					outfile.write('0')
					#sys.stdout.write( '0')
				if (bin2 != max_bin[chr2]):
					outfile.write('\t')
			outfile.write('\t')
		outfile.write('\n')
outfile.close()