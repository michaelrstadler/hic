'''
Scans the diagonal of a HI-C matrix
and counts the linkages between that bin and the bins to the left and the right, with 
both the distance included, and the directly adjacent bins to skip, submitted by the user.
So for 1 kb bins, inputs of -w 10 -d 3 would add up all the linkages between the bin and
bins between 3 and 10 kb to its left and to its right. The L/R score is just the simple
log ratio of right to left. 

Added a weight for the total reads used to compute directionality for a given bin. The L/R
directionality score is scaled by a, such that weighted directionality = (raw directionality)^a.

Input file is the "reduced diagonal" representation: chr	Lmost	Rmost	count

Writes a WIG file


Modified to allow it to also output a simple count of hte total Hi-C reads in the window (no direction component)
by using the -t --total tag
'''
from optparse import OptionParser
from math import log
import sys
import re
import gzip

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-w", "--width", dest="width",
					  help="width about diagonal, in bins", metavar="WIDTH")
	parser.add_option("-d", "--diagonal_skip", dest="diagonal_skip",
					  help="width about diagonal to skip, in bins", metavar="SKIP")
	parser.add_option("-a", "--scale_factor", dest="scale_factor", default=1,
					  help="Factor by which to scale weighting based on read count, w = [read_total]^a", metavar="SCALE")
	parser.add_option("-t", "--total", dest="total", default=False,
					  help="total mode, prints just the total reads in area defined by w and d, no direction", metavar="TOTAL")
	parser.add_option("-s", "--outstem", dest="outstem",
					  help="outfile stem", metavar="OUTSTEM")

	(options, args) = parser.parse_args()
	return options

# Initializes container based on the fly chromosome sizes (Dm3)            
def MakeBins(chromosomes, sizes):
	bin_counts = {}
	bin_size = int(options.bin_size)
	width = int(options.width)
	for i in range(0,6):
		num_bins = int(sizes[i] / bin_size)
		chr = chromosomes[i]
		bin_counts[chr] = {}
		for j in range(0, num_bins + 1):
			bin_counts[chr][j] = {}
	return bin_counts

# Reads reduced format into bin container
def Read_map(filename, chromosomes, sizes):
	bin_counts = MakeBins(chromosomes, sizes)
	filestring = str(options.filename)
	files = filestring.split(',')
	for file in files:
		print(file)
		if (file[-2:] == 'gz'):
			infile = gzip.open(file, 'rt')
		else:
			infile = open(file, 'r')
		for line in infile:
			line = line.rstrip()
			if (line[0] != '#'): #column/row totals ignore
				(chr, bin1, bin2, count) = line.split('\t')
				bin1 = int(bin1)
				bin2 = int(bin2)
				if (count == "NA"): count = 0
				if (abs(bin1 - bin2) < (2*int(options.width)) + 2): #just storing an extra couple bins to avoid thinking about end problems
					bin_counts[chr][bin1][bin2] = float(count)
		infile.close()
	return bin_counts

# Runs through and scores each genomic position by it's left/right deal
def Score_boundaries(bin_counts, sizes, chromosomes, scale_factor):
	bin_size = int(options.bin_size)
	width = int(options.width)
	diagonal_skip = int(options.diagonal_skip)
	outfilename = re.sub('.txt', '', options.outstem) + '_weightedDirectionality_w' + str(width) + 'd' + str(diagonal_skip) + 'a' + str(scale_factor) + '.WIG'
	outfile = open(outfilename,'w')
	name = 'width' + str(width) + '_skip' + str(diagonal_skip) + '_sf' + str(scale_factor)
	outfile.write ('track type=wiggle_0 name="' + name + '" description="test"' + '\n')
	for i in range(0,6):
		chr = chromosomes[i]
		max_bin = int(sizes[i] / bin_size)
		for bin1 in range(0, max_bin + 1):
			left_sum = 0
			right_sum = 0
			for i in range(bin1 - width, bin1 - diagonal_skip):
				if(i in bin_counts[chr][bin1]):
					left_sum += bin_counts[chr][bin1][i]
			for j in range(bin1 + diagonal_skip + 1, bin1 + width + 1):
				if(j in bin_counts[chr][bin1]):
					right_sum += bin_counts[chr][bin1][j]
			diff = log((right_sum + 0.5) / (left_sum + 0.5), 10) #pseudocounts
			weight = (left_sum + right_sum) ** scale_factor
			weighted_diff = diff * weight
			left_coord = str(bin1 * bin_size)
			right_coord = str((bin1 * bin_size) + bin_size - 1)
			if (options.total):
				outfile.write('chr' + chr + '\t' + left_coord + '\t' + right_coord + '\t' + str(left_sum + right_sum) + '\n')
			else:
				outfile.write('chr' + chr + '\t' + left_coord + '\t' + right_coord + '\t' + str(weighted_diff) + '\n')
	outfile.close()

# Main section calling everything
#dm3
#sizes = (23011544, 22422827, 24543557, 1351857, 21146708, 27905053)
#dm6
sizes = (23513712, 23542271, 28110227,1348131,25286936,32079331)
chromosomes = ('2L','X','3L','4','2R','3R')

options = parse_options()
print("Reading map from:")
bin_counts = Read_map(options.filename, chromosomes, sizes)
print("Done reading map...")
scale_factor = float(options.scale_factor)
Score_boundaries(bin_counts, sizes, chromosomes, scale_factor)
print("Done calculating directionality!")

