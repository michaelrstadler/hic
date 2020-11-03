"""

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
					  
	parser.add_option("-s", "--species",
						dest="species", default='mel',
					  help="species: 'mel' 'pse'")
	
	parser.add_option("-o", "--output_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds _binCounts_Xkb.txt")

	(options, args) = parser.parse_args()
	return options

def get_chromosome_max_bins(species, bin_size):
	chromosome_bin_numbers = {}
	if (species == 'mel'):
			chromosome_bin_numbers = {
				"2L":int(23513712 / bin_size),
				"2R":int(25286936 / bin_size),
				"3L":int(28110227 / bin_size),
				"3R":int(32079331 / bin_size),
				"4":int(1348131 / bin_size),
				"Y":int(3667352 / bin_size),
				"X":int(23542271 / bin_size),
			}

	if (species == 'pse'):
		chromosome_bin_numbers = {
			"2":int(30711475 / bin_size),
			"3":int(19738957 / bin_size),
			"2":int(30711475 / bin_size),
			"2":int(30711475 / bin_size),
			"2":int(30711475 / bin_size),
			"2":int(30711475 / bin_size),
		}	
	return(chromosome_bin_numbers)
                
def make_bin_container(chromosome_max_bins):
	bin_container = {}

	for chr1 in chromosome_max_bins.keys():
		bin_container[chr1] = {}
		#print(chr1)
		for bin1 in range(0, chromosome_max_bins[chr1] + 1):
			bin_container[chr1][bin1] = {}
			for chr2 in chromosome_max_bins.keys():
				bin_container[chr1][bin1][chr2] = {}
				for bin2 in range(0, chromosome_max_bins[chr2] + 1):
					bin_container[chr1][bin1][chr2][bin2] = 0
	return(bin_container)



options = parse_options()
bin_size = int(options.bin_size)
species = options.species
max_bin = get_chromosome_max_bins(species, bin_size)
bin_bin_counts = make_bin_container(max_bin)
#print(bin_bin_counts['2L'][1]['3R'][10])
#quit()
filenames = options.filenames
files = filenames.split(',')
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
		bin1 = int(Lmost1 / bin_size)
		bin2 = int(Lmost2 / bin_size) 
		if (chr1 in max_bin and chr2 in max_bin):
			#print(chr1 + ' ' + str(bin1) + ' ' + chr2 + ' ' + str(bin2))
			bin_bin_counts[chr1][bin1][chr2][bin2] = bin_bin_counts[chr1][bin1][chr2][bin2] + 1
			# Avoid double counting on diagonal. This will be triggered unless chromosome and bin are the same
			if (chr1 != chr1 or bin1 != bin2):
				bin_bin_counts[chr2][bin2][chr1][bin1] = bin_bin_counts[chr2][bin2][chr1][bin1] + 1

	file1.close()

file_stem = ''

if (species == 'mel'):
	chromosomes = ('X', '2L', '2R', '3L', '3R', '4');

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

outfile = open(file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt','w')

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
				outfile.write(str(count))
				if (bin2 != max_bin[chr2]):
					outfile.write('\t')
			outfile.write('\t')
		outfile.write('\n')
outfile.close()