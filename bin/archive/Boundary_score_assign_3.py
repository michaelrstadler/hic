'''
Takes a file of directionality scores, computes the average directionality within some window to the left and right
of each bin, reports these averages as a WIG-style file. This is a bit of a legacy as it once produced a true WIG file.
The output of this script can be read by R script boundaries.call to call boundaries with different thresholds. 
'''
from optparse import OptionParser
from math import log
import numpy as np
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Directional bias file", metavar="FILE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-w", "--width", dest="width",
					  help="width about diagonal, in bins", metavar="WIDTH")
	parser.add_option("-z", "--zeroes_tolerated", dest="zeros",
					  help="number of zeros allowed in window to have a score", metavar="ZEROS")
	parser.add_option("-n", "--track_name", dest="track_name",
					  help="name for browser track", metavar="TRACKNAME")

	(options, args) = parser.parse_args()
	return options

# Adds up counts within window specified
def sum_counts(chr, start, end, counts, bin_size):
	sum = 0
	for i in range(start, end + 1):
		bin = i * bin_size
		if (bin in counts[chr]):
			sum += counts[chr][bin]
	return(sum)

def update_max(max_bins, chr, bin_number):
	if (chr in max_bins):
		if(bin_number > max_bins[chr]):
			max_bins[chr] = bin_number
	else:
		max_bins[chr] = bin_number

def add_count(chr, bin, count, counts):
	if (chr in counts):
		counts[chr][bin] = count
	else:
		counts[chr] = {}
		counts[chr][bin] = count

def count_zeros(chr, bin, width, counts, bin_size):
	zero_count = 0
	for i in range(bin - width, bin + width - 1):
		bin_local = i * bin_size
		if (bin_local in counts[chr]):
			if (counts[chr][bin_local] == 0.0):
				zero_count +=1
	return zero_count


options = parse_options()
bin_size = int(options.bin_size)
max_zeros = int(options.zeros)
counts = {}
max_bins = {}
infile = open(options.filename,'r')
for line in infile:
	line = line.rstrip()
	if (line[0:5] != 'track'):
		items = line.split('\t')
		(chr, bin, end, count) = items
		bin = int(bin)
		count = float(count)
		bin_number = bin / bin_size # bit of confusion here betwee "bin number", which is 0, 1, 2, and the positions of the bins, 0, 500, 1000...
		update_max(max_bins, chr, bin_number)
		add_count(chr, bin, count, counts)
		
outfile_stem = re.sub('.WIG', '', options.filename)
outfile = open(outfile_stem + '_boundarieScore_w' + options.width  + '.WIG','w')
#outfile.write('track type=wiggle_0 name="boundarieScore_' + options.track_name + '" description="boundarieScore_' + options.track_name + '"\n')
chromosomes = ('2L','X','3L','4','2R','3R')
scores = []

for chr in chromosomes:
	chr = 'chr' + chr
	max_bin = int(max_bins[chr])
	width = int(options.width)
	for bin1 in range(0, max_bin + 1):
		left_sum = sum_counts(chr, bin1 - width, bin1 - 1, counts, bin_size) 
		right_sum = sum_counts(chr, bin1, bin1 + width - 1, counts, bin_size)
		bin_pos = bin1 * bin_size
		score = 0
		if (count_zeros(chr, bin1, width, counts, bin_size) <= max_zeros):
			left_average = left_sum / width
			right_average = right_sum / width
			diff = right_average - left_average
			scores.append(chr + '\t' + str(bin1) +  '\t' + str(diff))


for line in scores:
	outfile.write(line + '\n')
		

outfile.close()
infile.close()

