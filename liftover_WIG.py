'''
Implements PyLiftover to a 4 column BED file. The problem that this solves specifically is that in lifting
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
					  help="new genome", metavar="NEW")
	(options, args) = parser.parse_args()
	return options

options = parse_options()
chrom_max_size = 4e7
bin_size = 20

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
		items = line.split()
		if (len(items) >= 3):
			[chr1, pos1, pos2, val] = [items[i] for i in [0,1,2,3]]
			try:
				pos1 = int(pos1)
				pos2 = int(pos2)
			except:
				outfile.write(line + '\n')
				pass
			if (chr1[0:3] != 'chr'):
				chr1 = 'chr' + chr1
			liftover1 = lo.convert_coordinate(chr1, pos1)
			liftover2 = lo.convert_coordinate(chr1, pos2)
			items_liftedover = items
			if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
				if (isinstance(liftover2, list) and len(liftover2) > 0 and len(liftover2[0]) >= 2):
					liftover_bin1 = liftover1[0][1]
					liftover_bin2 = liftover2[0][1]

					for pos in range(liftover_bin1, liftover_bin2 + 1):
						bin = int(pos / bin_size)
						if (chr1 in chrom_max_bin):
							if (bin > chrom_max_bin[chr1]):
								chrom_max_bin[chr1] = bin
							bin_counts[chr1][bin] = bin_counts[chr1][bin] + float(val)
						else:
							chrom_max_bin[chr1] = bin
							bin_counts[chr1] = zeros(int (chrom_max_size / bin_size))
		else:
			outfile.write(line + '\n')

	for chr2 in chrom_max_bin:
		in_zero_run = False
		zero_start_bin = 0
		for bin2 in range (0, chrom_max_bin[chr2] + 1):
			pos3 = bin2 * bin_size
			pos4 = pos3 + bin_size - 1
			value = bin_counts[chr2][bin2] / bin_size
			if (value > 0):
				# check to see if we have trailing zeros
				if (in_zero_run):
					outfile.write(chr2 + '\t' + str(zero_start_bin * bin_size) + '\t' + str(pos3 - 1) + '\t0\n')
				#print this line
				outfile.write(chr2 + '\t' + str(pos3) + '\t' + str(pos4) + '\t' + str(value) + '\n')
				in_zero_run = False
			else:
				if(not in_zero_run):
					in_zero_run = True
					zero_start_bin = bin2
	infile.close()
	outfile.close()
