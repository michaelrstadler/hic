# -*- coding: utf-8 -*-
"""
Takes two bowtie outputs from single end alignments of paired reads, matches
based on read name and prints paired end alignments in format:

read name \t R1_strand \t R1_chromosome \t R1_5'pos \t R2_strand \t R2_chromosome \t R2_5'pos

Change in version 3 is that it reports the position of the 5'-most bp in the read. This compensates
for decreases in quality at the 3' that result in mapping variability and artifacts when removing
PCR duplicates.

@author: MStadler
"""
import sys
from optparse import OptionParser
import re
import gzip

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file1", dest="filename1",
					  help="First mapping file", metavar="FILE1")
	parser.add_option("-s", "--file2", dest="filename2",
					  help="Second mapping file", metavar="FILE2")

	(options, args) = parser.parse_args()
	return options

def Get5prime(splitline):
	strand = splitline[1]
	chr = splitline[2]
	pos = splitline[3]
	seq = splitline[4]
	if (strand == '+'):
		pass
	else:
		size = len(seq)
		pos = str(int(pos) + size)
	return(pos)
		
options = parse_options()

file1_data = {}
goodPair_count = 0
total = 0
distance_minimum = 0
outfilename = ''
if (re.search('R1', options.filename1)):
	outfilename = options.filename1
elif (re.search('R1', options.filename2)):
	outfilename = options.filename2
outfilename = re.sub('_R1_','_', outfilename)
outfilename = re.sub('.gz', '', outfilename)
outfilename = re.sub('.bowtie', '', outfilename)
outfilename = outfilename +  '_pairMergedV3.txt'
outfile = open(outfilename, "w")

f1 = options.filename1
if (f1[-2:] == 'gz'):
		file1 = gzip.open(f1, 'rt')
else:
	file1 = open(f1, 'r')

for line in file1:
    pair_name = line.split()[0]
    splitline= line.split('\t')
    pos5p = Get5prime(splitline)
    data = "\t".join(splitline[1:3] + [pos5p])
    file1_data[pair_name] = data
file1.close()

f2 = options.filename2
if (f2[-2:] == 'gz'):
		file2 = gzip.open(f2, 'rt')
else:
	file2 = open(f2, 'r')

for line in file2:
	total += 1
	pair_name = line.split()[0]
	splitline2= line.split('\t')
	if pair_name in file1_data:
		splitline1 = file1_data[pair_name].split('\t')
		pos_5p_2 = Get5prime(splitline2)
		outfile.write(pair_name + "\t" + splitline1[0] + '\t' + splitline1[1] + '\t' + splitline1[2] + '\t' + splitline2[1] + '\t' + splitline2[2] + '\t' + pos_5p_2 + '\n')
		goodPair_count += 1

file2.close()


pct = float(goodPair_count) / float(total) * 100
#print >> sys.stderr, str(goodPair_count) + ' of ' + str(total) + ' reads have aligned pairs, or ' + str(pct) + ' percent'
