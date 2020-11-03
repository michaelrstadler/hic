'''
Implements PyLiftover to a file in which the first three columns of every line are chr, pos1, pos2, all tab-delimited. 
Just replaces the two positions with lifted-over positions and leave the rest intact.

If the first three columns don't seem correct (string, int, int), just prints the line.
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
	(options, args) = parser.parse_args()
	return options

options = parse_options()

lo = LiftOver(options.original, options.new)

f1 = options.infile
if (f1[-2:] == 'gz'):
	infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

for line in infile:
	line = line.rstrip()
	items = line.split()
	if (len(items) >= 3):
		[chr1, pos1, pos2] = [items[i] for i in [0,1,2]]
		try:
			pos1 = int(pos1)
			pos2 = int(pos2)
		except:
			print(line)
			pass
		if (chr1[0:3] != 'chr'):
			chr1 = 'chr' + chr1
		liftover1 = lo.convert_coordinate(chr1, pos1)
		liftover2 = lo.convert_coordinate(chr1, pos2)
		items_liftedover = items
		if (isinstance(liftover1, list) and len(liftover1) > 0 and len(liftover1[0]) >= 2):
			if (isinstance(liftover2, list) and len(liftover2) > 0 and len(liftover2[0]) >= 2):
				#items_liftedover = items[0] + list(liftover1[0][1]) + list(liftover2[0][1]) + items[3:]  #items[0:2] + list(liftover1[0][0:2]) + list(items[4]) + list(liftover2[0][0:2])
				#items_liftedover = list(map(str, items_liftedover))
				#print('\t'.join(items_liftedover))
				print(items[0] + '\t' + str(liftover1[0][1]) + '\t' + str(liftover2[0][1]) + '\t' + '\t'.join(items[3:]))
	else:
		print(line)

infile.close()