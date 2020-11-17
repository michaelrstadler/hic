'''
Converts a BED format file to a gff3-like file for input into scripts that require it.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--bed_file", dest="bed_file",
					  help="BED file", metavar="BEDFILE")
	(options, args) = parser.parse_args()
	return options

options = parse_options()
infile = open(options.bed_file, 'r')
outfile = open(re.sub('bed','gff3', options.bed_file),'w')

for line in infile:
	line = line.rstrip()
	items = line.split()
	outfile.write(items[0] + '\t.\t' + items[3] + '\t' + items[1] + '\t' + items[2] + '\t' + items[4] + '\t.\t.\t.\n')

infile.close()
outfile.close()