'''
Implements PyLiftover to a file tab-delimited file in which some columns represent the chromosome, left position and right position. 
User indicates columns (0-based) that contain the relevant features.
Just replaces the two positions with lifted-over positions and leave the rest intact.
'''
from optparse import OptionParser
import sys
import re
from pyliftover import LiftOver
import gzip


def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="infile",
					  help="input file", metavar="INFILE")
	parser.add_option("-o", "--original", dest="original",
					  help="original genome", metavar="ORIGINAL")
	parser.add_option("-n", "--new", dest="new",
					  help="new genome", metavar="NEW")
	parser.add_option("-c", "--chr", dest="chr",
					  help="column with chromosome (0-based)", metavar="CHR")
	parser.add_option("-l", "--left", dest="left",
					  help="column with left position (0-based)", metavar="LEFT")
	parser.add_option("-r", "--right", dest="right",
					  help="column with right position (0-based)", metavar="RIGHT")
	(options, args) = parser.parse_args()
	return options

options = parse_options()

lo = LiftOver(options.original, options.new)

f1 = options.infile
chr_col = int(options.chr)
left_pos_col = int(options.left)
right_pos_col = int(options.right)

if (f1[-2:] == 'gz'):
	infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

for line in infile:
	items = line.split()
	pos1 = int(items[left_pos_col])
	pos2 = int(items[right_pos_col])
	chr = items[chr_col]
	if (chr[0:3] != 'chr'):
		chr = 'chr' + chr
	liftover1 = lo.convert_coordinate(chr, pos1)
	liftover2 = lo.convert_coordinate(chr, pos2)
	items_liftedover = items
	items_liftedover[left_pos_col] = str(liftover1[0][1])
	items_liftedover[right_pos_col] = str(liftover2[0][1])
	#complicated error checking to make sure liftover worked. Can't remember why this is such a mess but it is.
	if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
		if (isinstance(liftover2, list) and len(liftover2) > 0 and len(liftover2[0]) >= 2):
			print('\t'.join(items_liftedover))
		else:
			print("Did not liftover")
	else:
		print('Did not liftover')

infile.close()