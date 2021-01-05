'''
Implements PyLiftover to a Hi-C viewer track file. The problem that this solves specifically is that in lifting
over, the resulting lines are often not in order and even occasionally overlap. This script solves that by
reading the data in for every nucleotide position in the new coordinates and reprinting the data for each 
chromosome in order. In the case of overlapping values, it arbitrarily takes whatever is latest in the file.

To avoid slow step of re-loading liftover data, written to handle an entire folder of files to lift.
'''
from optparse import OptionParser
import sys
import re
from pyliftover import LiftOver
import gzip
from os import listdir
from os.path import isfile, join
from numpy import zeros


def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--folder", dest="infolder",
					  help="folder containing input files", metavar="INFOLDER")
	parser.add_option("-o", "--original", dest="original",
					  help="original genome", metavar="ORIGINAL")
	parser.add_option("-n", "--new", dest="new",
					  help="new genome", metavar="BINSIZE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size", default=500)

	(options, args) = parser.parse_args()
	return options

options = parse_options()
chrom_max_size = 4e7
bin_size = int(options.bin_size)

lo = LiftOver(options.original, options.new)

for f in listdir(options.infolder):
	print(f)
	if (f[0] == '.'): #takes care of .ds_store file
		continue
	f_path = join(options.infolder, f)
	bin_counts = {}
	chrom_max_bin = {}

	if (f_path[-2:] == 'gz'):
		infile = gzip.open(f_path, 'rt')
	else:
		infile = open(f_path, 'r')
	
	out_stem = re.sub('.gz|.txt|.wig|.WIG', '', f)
	outfile = open(join(options.infolder, out_stem) + '_' + options.new + '.txt', 'w')

	for line in infile:
		line = line.rstrip()
		[chr1, bin_, val] = line.split()
		
		if (chr1[0:3] != 'chr'):
			chr1 = 'chr' + chr1

		pos1 = int(bin_) * bin_size
		liftover1 = lo.convert_coordinate(chr1, pos1)
		if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
			
			pos_new = liftover1[0][1]
			bin_new = int(pos_new // bin_size)
			outfile.write(chr1 + '\t' + str(bin_new) + '\t' + val + '\n')
					
	infile.close()
	outfile.close()
