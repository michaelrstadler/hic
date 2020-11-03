'''
Takes "diagonal bin coverage" from HiC_count_bin_linkages_singleChr_diagonal_v1a.py. This
is a compressed way to store information that only stores information within some radius of
the diagonal. It also takes a list of positions (tab-delimited: name\tchr\tposition) and
a width and generate local matrices +- the width selected of each position.

v1b: more efficiently stores data (uses way less memory)


'''
from optparse import OptionParser
import sys
import re
import gzip

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-l", "--locations", dest="loc_file",
					  help="File of locations to map", metavar="FILE2")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-s", "--file_stem", dest="file_stem",
					  help="stem for outfiles", metavar="STEM")
	parser.add_option("-w", "--width", dest="width",
					  help="num bins on either side of diagonal to store, in bins", metavar="WIDTH")

	(options, args) = parser.parse_args()
	return options
            
def MakeBins():
	bin_counts = {}
	chromosomes = ('2L','X','3L','4','2R','3R')
	bin_size = int(options.bin_size)
	sizes = (23011544, 22422827, 24543557, 1351857, 21146708, 27905053)
	for i in range(0,6):
		num_bins = int(sizes[i] / bin_size)
		chr = chromosomes[i]
		bin_counts[chr] = {}
		for j in range(0, num_bins + 1):
			bin_counts[chr][j] = {}
	return bin_counts

def Generate_map(name, chr, pos1, pos2):
	outfile_local = options.file_stem + '_' + name + '.txt'
	outfile = open(outfile_local, 'w')
	#print colnames
	for k in range (pos1, pos2 + 1):
		if (k != pos1):
			outfile.write('\t')
		outfile.write(chr + '_' + str(k))
	outfile.write('\n')
		
	for i in range (pos1, pos2 + 1):
		outfile.write(chr + '_' + str(i)) #row name
		for j in range(pos1, pos2 + 1):
			count = '0.0'
			if (i in bin_counts[chr]):
				if(j in bin_counts[chr][i]):
					count = bin_counts[chr][i][j]
			outfile.write('\t' + count)
		outfile.write('\n')
			
	outfile.close()			

def add_count(bin_counts, chr, bin1, bin2, count):
	if (chr in bin_counts):
		if (bin1 in bin_counts[chr]):
			bin_counts[chr][bin1][bin2] = count
		else:
			bin_counts[chr][bin1] = {}
			bin_counts[chr][bin1][bin2] = count
	else:
		bin_counts[chr] = {}
		bin_counts[chr][bin1] = {}
		bin_counts[chr][bin1][bin2] = count

options = parse_options()
#bin_counts = MakeBins()
bin_counts = {}
width_to_store = int(options.width)

f1 = options.filename
if (f1[-2:] == 'gz'):
		infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

for line in infile:
	line = line.rstrip()
	if (line[0] != '#'):
		(chr, bin1, bin2, count) = line.split('\t')
		bin1 = int(bin1)
		bin2 = int(bin2)
		if (abs(bin1 - bin2) < (width_to_store + 2)): #just storing an extra couple bins to avoid thinking about end problems
			add_count(bin_counts, chr, bin1, bin2, count)
		#bin_counts[chr][bin1][bin2] = count

infile.close()

locfile = open(options.loc_file, 'r')

for line in locfile:
	line = line.rstrip()
	items = line.split('\t')
	name = items[0]
	chr = items[1]
	chr = re.sub('chr', '', chr)
	pos1 = items[2]
	pos2 = items[3]
	bin1 = int(int(pos1) / int(options.bin_size))
	bin2 = int(int(pos2) / int(options.bin_size))
	#if ((bin1 + 1) in bin_counts[chr][bin1]): #this is just a test for whether this chromosome was loaded. Allows inclusion of all locations in single file
	#	Generate_map(name, chr, bin1, bin2)
	if (chr in bin_counts):
		Generate_map(name, chr, bin1, bin2)

locfile.close()

