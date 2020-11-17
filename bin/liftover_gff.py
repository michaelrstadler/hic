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
bin_size = 10

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
	
	out_stem = re.sub('.gz|.txt|.wig|.WIG|.gff3', '', f)
	outfile = open(join(options.infolder, out_stem) + '_' + options.new + '.gff3', 'w')

	for line in infile:
		line = line.rstrip()
		if (line[0] != '#'):
			items = line.split()
			[chr1, pos1, pos2] = [items[i] for i in [0,3,4]]
			try:
				pos1 = int(pos1)
				pos2 = int(pos2)
			except:
				outfile.write(line + '\n')
				pass
			if (options.original == 'dm3' and chr1[0:3] != 'chr'):
				chr1 = 'chr' + chr1
			
			liftover1 = lo.convert_coordinate(chr1, pos1)
			liftover2 = lo.convert_coordinate(chr1, pos2)
			items_liftedover = items
			if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
				if (isinstance(liftover2, list) and len(liftover2) > 0 and len(liftover2[0]) >= 2):
					liftover_bin1 = liftover1[0][1]
					liftover_bin2 = liftover2[0][1]

					items_liftedover[3] = str(liftover_bin1)
					items_liftedover[4] = str(liftover_bin2)
					line_liftedover = ('\t').join(items_liftedover)
					outfile.write(line_liftedover + '\n')
		else:
			outfile.write(line + '\n')

	infile.close()
	outfile.close()
