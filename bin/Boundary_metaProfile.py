'''
Takes a file of positions (e.g., boundaries) in UCSC format and a wig file of some genomic feature (ChIP), sums
the wig feature relative to the left-most position of the boundaries, reports the aggregate profile. I built this
normalization function, which divides the row generated around each boundary by the sum...but it works poorly 
and mathematically doesn't make a lot of sense. Don't feel like explaining, but I left it just in case it's
useful somehow. Anyway don't use the -n function.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-p", "--positions_file", dest="pos_file",
					  help="file with positions", metavar="POSFILE")
	parser.add_option("-f", "--features_file", dest="features_file",
					  help="files with features, comma-delimited", metavar="FEATFILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile", metavar="OUTFILE")
	parser.add_option("-w", "--window", dest="window", default=20,
					  help="window in bins", metavar="WINDOW")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100,
					  help="bin size bp", metavar="BINSIZE")
	parser.add_option("-i", "--individual", dest="individual", default=False,
					  help="Print individual lines? (any character)", metavar="INDIV")
	parser.add_option("-g", "--gff", dest="gff", default=False,
					  help="Is position file in gff format?", metavar="GFF")

	(options, args) = parser.parse_args()
	return options

def add_position(line, genome, window, bin_size, counts, individual, outfile):
	line = line.rstrip()
	(chr, pos1, pos2) = ('','','')
	if (options.gff):
		items = line.split()
		if (len(items) >= 5):
			chr = items[0]
			pos1 = items[3]
			pos2 = items[4]
	else:
		chr, pos1, pos2 = line.split('\t')[0:3]
		#(pos1, pos2) = positions.split('-')
	bin_Lmost = int(int(pos1) / bin_size)
	local_counts = {}
	local_sum = 0

	if (not re.match('chr', chr)): chr = 'chr' + chr

	if (chr in genome):
		for i in range(-1 * window, window):
			bin = bin_Lmost + i
			if (bin in genome[chr]):
				local_counts[i] = genome[chr][bin]
				local_sum += genome[chr][bin]
			else:
				local_counts[i] = 0
		if (individual):
			outfile.write(chr + ':' + pos1 + '-' + pos2)
			for i in range(-1 * window, window):
				outfile.write('\t' + str(local_counts[i]))
			outfile.write('\n')
		add_row(local_counts, counts, window)

def normalize_row(local_counts, local_sum):
	norm_counts = {}
	for i in local_counts:
		norm_counts[i] = local_counts[i] / local_sum
	return(norm_counts)

def add_row(local_counts, counts, window):
	for i in range(-1 * window, window):
		if (i in local_counts):
			counts[i] += local_counts[i]
	

def add_feature(line, genome, bin_size):
	(chr, pos1, pos2, value) = line.split()
	pos1 = int(pos1)
	pos2 = int(pos2)
	value = float(value)
	bin = int(pos1 / bin_size)

	if (chr in genome):
		if (bin in genome[chr]):
			genome[chr][bin] += value
		else:
			genome[chr][bin] = value
	else:
		genome[chr] = {}
		genome[chr][bin] = value
			


options = parse_options()
outfile = open(options.outfile,'w')
genome = {}
counts = {}
window = int(options.window)
bin_size = int(options.bin_size)
individual = options.individual
for i in range(-1 * window, window):
	counts[i] = 0

# first load feature file WIG
featurefile = open(options.features_file, 'r')

for line in featurefile:
	line = line.rstrip()
	if (re.search('track', line)):
		pass
	else:
		add_feature(line, genome, bin_size)

featurefile.close()

# Next deal with positions file
pos_file = open(options.pos_file, 'r')

for line in pos_file:
	if (line[0] != '#'):
		add_position(line, genome, window, bin_size, counts, individual, outfile)

pos_file.close()

if (not individual):
	for i in range(-1 * window, window):
		outfile.write(str(i) + '\t' + str(counts[i]) + '\n')

outfile.close()




