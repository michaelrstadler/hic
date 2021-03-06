"""
 Takes BED files of genomic features (such as promoters) and generates a 4-column WIG-style
 file for the Hi-C viewer, with 1 where features are and 0 for other positions.

 Wrote to parse files from promoter database https://epd.epfl.ch/

"""

from optparse import OptionParser
import sys
import re
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="bed file", metavar="FILE")
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
	parser.add_option("-n", "--track_name",
					   dest="track_name", default='temp_name',
					  help="track name")
	parser.add_option("-d", "--track_desc",
					   dest="track_desc", default='temp_desc',
					  help="track description")
	parser.add_option("-o", "--output",
					   dest="output", default='v',
					  help="output format: w for wig, v for viewer track")


	(options, args) = parser.parse_args()
	return options


def print_wig(bin_counts, bin_size, track_name, track_desc, max_bin, outfile):
	outfile.write('track type=wiggle_0 name="' + track_name + '" description="' + track_desc + '"\n')
	for chrom in bin_counts:
		outfile.write('variableStep chrom=chr' + chrom + ' span=' + str(bin_size) + '\n')
		for pos in range(0, max_bin[chrom]):
			if (pos in bin_counts[chrom]):
				outfile.write(str((pos * bin_size) + 1) + '\t' + str(bin_counts[chrom][pos]) + '\n')
			else:
				outfile.write(str((pos * bin_size) + 1) + '\t' + '0' + '\n')

def print_viewer_track(bin_counts, bin_size, max_bin, outfile):
	for chrom in bin_counts:
		for pos in range(0, max_bin[chrom]):
			start = str(pos * bin_size)
			end = str(int(start) + bin_size - 1)
			if (pos in bin_counts[chrom]):
				outfile.write(chrom + '\t' + start + '\t' + end + '\t' + str(bin_counts[chrom][pos]) + '\n')
			else:
				outfile.write(chrom + '\t' + start + '\t' + end + '\t' + '0' + '\n')

options = parse_options()
bin_size = int(options.bin_size)
max_bin = {}
counts = {}
track_name = options.track_name
track_desc = options.track_desc

file1 = open(options.filename, 'r')

for line in file1:
	line = line.rstrip()
	if (line[0:3] == 'chr'):
		items = line.split('\t')
		chrom = items[0] 
		chrom = re.sub('chr', '', chrom) # takes off leading 'chrom'
		pos = int(items[1])
		binnum = int(int(pos) / bin_size)
		value = 1
		if (items[5] == '-'):
			value = -1 #placeholder for BED...can modulate this if desired later
		if(chrom in counts):
			if (binnum in counts):
				counts[chrom][binnum] = counts[chrom][binnum] + value
			else:
				counts[chrom][binnum] = value
		else:
			counts[chrom] = {}
			counts[chrom][binnum] = value
		
		if(chrom in max_bin):
			if (binnum > max_bin[chrom]):
				max_bin[chrom] = binnum
		else:
			max_bin[chrom] = binnum

file1.close()
file_stem = re.sub('.txt', '', options.filename)

chromomosomes = ['X','2L','2R','3L','3R','4']

if (options.output == 'w'):
	outfile = open(file_stem + '_binCounts_' + str(bin_size) + 'bp.wig','w')
	print_wig(counts, bin_size, track_name, track_desc, max_bin, outfile)	

elif (options.output == 'v'):
	outfile = open(file_stem + '_binCounts_' + str(bin_size) + 'bp.txt','w')
	print_viewer_track(counts, bin_size, max_bin, outfile)



	
outfile.close()



