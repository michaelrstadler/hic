'''
Implements PyLiftover to convert a paired mapping file from dm6 coordinates to dm3 coordinates.
'''
from optparse import OptionParser
import sys
import re
from pyliftover import LiftOver
import gzip


def parse_options():
	usage = 'Converts pairmap file from dm6 coordinates to dm3 coordinates'
	parser = OptionParser(usage=usage)
	parser.add_option("-p", "--pairfile", dest="infile",
					  help="input file: paired mapped", metavar="INFILE")
	(options, args) = parser.parse_args()
	return options


options = parse_options()

lo = LiftOver('dm6', 'dm3')

f1 = options.infile
if (f1[-2:] == 'gz'):
	infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

for line in infile:
	items = line.split()
	[chr1, pos1, chr2, pos2] = [items[i] for i in [2,3,5,6]]
	pos1 = int(pos1)
	pos2 = int(pos2)
	chr1 = 'chr' + chr1
	chr2 = 'chr' + chr2
	liftover1 = lo.convert_coordinate(chr1, pos1)
	liftover2 = lo.convert_coordinate(chr2, pos2)
	items_liftedover = items
	if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
		if (isinstance(liftover2, list) and len(liftover2) > 0 and len(liftover2[0]) >= 2):
			items_liftedover = items[0:2] + list(liftover1[0][0:2]) + list(items[4]) + list(liftover2[0][0:2])
			items_liftedover = list(map(str, items_liftedover))
			print('\t'.join(items_liftedover))

infile.close()