'''
Cleans up WIG files for use in my scripts that require 4-column WIG files with single track lines. Keeps only the first
track line, converts a two-column (fixed span) WIG file to a 4-column WIG format, leaves 4-column intact
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
outfilename = ''
if re.search('WIG$',options.wig_file):
	outfilename = re.sub('.WIG$','_4col.WIG', options.wig_file)
elif re.search('wig$',options.wig_file):
	outfilename = re.sub('.wig$','_4col.WIG', options.wig_file)
else:
	sys.exit('Does not appear to be a WIG file!')
outfile = open(outfilename, 'w')


curr_chr = ''
curr_span = 0
curr_pos = 0
track_line_written = False


for line in infile:
	line = line.rstrip()
	if (line[0:5] == 'track'):
		if(not track_line_written):
			outfile.write(line + '\n')
			track_line_written = True
	elif(line[0:4] == 'fixe'):
		items = line.split()
		curr_chr = items[1]
		curr_chr = curr_chr[6:]
		curr_pos = 0
		curr_span = items[4]
		curr_span = int(curr_span[5:])
		if (not re.match('chr', curr_chr)):
			curr_chr = 'chr' + curr_chr
	elif(line[0] == '#'):
		pass
	else:
		line_values = line.split()
		if (len(line_values) == 4):
			outfile.write(line + '\n')
		elif(len(line_values) == 2):
			(pos, value) = line.split()
			pos = int(pos)
			outfile.write(curr_chr + '\t' + str(pos) + '\t' + str(pos + curr_span - 1) + '\t' + value + '\n')
infile.close()