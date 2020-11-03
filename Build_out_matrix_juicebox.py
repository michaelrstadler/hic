'''
Takes Juicebox matrix "dump" file and builds it out into a proper matrix with column and row names in standard format
"chr_pos". The juicebox format is bin1 bin2 value. It does NOT repeat the reciprocal values (10 11 is there and 11 10 isn't)
so to make symmetric matrix, have to handle that. Also, there are a bunch of missing values. I don't know why, but it is
conceivable that they are 0s. I have the script printing the missing values and subbing them with 0, but I also tried
subbing with the average value (~1.3). Didn't seem to make a big difference.
'''
from optparse import OptionParser
import sys
import re
import random
from random import randint
from random import shuffle
from numpy import percentile

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="file",
					  help="JuiceBox file matrix dump", metavar="FILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile", metavar="FILE")
	parser.add_option("-c", "--chromosome", dest="chromosome",
					  help="chromosome ", metavar="CHR")
	parser.add_option("-b", "--bin_size", dest="bin_size",
					  help="bin size", metavar="BIN")
	parser.add_option("-n", "--num_bins", dest="num_bins",
					  help="num bins", metavar="NUMBINS")
	(options, args) = parser.parse_args()
	return options

def add_counts(counts, bin1, bin2, value):
	if (bin1 in counts):
		counts[bin1][bin2] = value
	else:
		counts[bin1] = {}
		counts[bin1][bin2] = value

def update_maxbin(max_bin, bin):
	if (bin > max_bin):
		return(bin)
	return(max_bin)


options = parse_options()
chr = options.chromosome
if (not re.match('chr', chr)):
	chr = 'chr' + chr
bin_size = int(options.bin_size)
counts = {}
max_bin = 0
num_bins = int(options.num_bins)

infile = open(options.file, 'r')
outfile = open(options.outfile, 'w')
for line in infile:
	line = line.rstrip()
	(bin1, bin2, value) = line.split()
	bin1 = int(int(bin1) / bin_size)
	bin2 = int(int(bin2) / bin_size)
	max_bin = update_maxbin(max_bin, bin1)
	max_bin = update_maxbin(max_bin, bin2)
	#print(max_bin)
	add_counts(counts, bin1, bin2, value)
	if (bin1 != bin2):
		add_counts(counts, bin2, bin1, value)

for i in range(0,num_bins):
	if (i != 0):
		outfile.write('\t')
	outfile.write(chr + '_' + str(i))
outfile.write('\n')

for i in range(0,num_bins):
	outfile.write(chr + '_' + str(i))
	for j in range(0, num_bins):
		if (i in counts and j in counts[i]): #seems to be a bug where Juicebox doesn't always dump every bin
			outfile.write('\t' + counts[i][j]) 
		else:
			#print(str(i) + '   ' + str(j))
			outfile.write('\t' + '0') #average is about 1.3, if I want to use that
	outfile.write('\n')


infile.close()
outfile.close()