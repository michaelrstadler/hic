'''
This takes modEncode WIG files and makes them UCSC viewable. First thing is it appends 'chr' to chromosome names
and also throws out all the "track" lines after the first one. This is important so taht all of the data is
in a single track instead of being useless split by chr.
'''
from optparse import OptionParser
import sys
import re
from Bio import SeqIO

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="wig_file",
					  help="WIG file", metavar="FILE")

	(options, args) = parser.parse_args()
	return options

options = parse_options()

infile = open(options.wig_file,'r')
outfilename = re.sub('.wig', '_ucsc.wig', options.wig_file)
outfile = open(outfilename, 'w')

curr_chr = ''
curr_span = 0
trackline_printed = False


for line in infile:
	line = line.rstrip()
	if (line[0:5] == 'track'):
		if(not trackline_printed):
			outfile.write(line + '\n')
			trackline_printed = True
		pass
	elif(line[0:4] == 'vari'):
		if (not re.search('chrom=chr', line)):
			print ('booya')
			line = re.sub('chrom=', 'chrom=chr',line)
		outfile.write(line + '\n')
	else:
		outfile.write(line + '\n')

infile.close()
outfile.close()