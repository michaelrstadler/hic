'''

'''
from optparse import OptionParser
import sys
import re
import gzip
import random
import subprocess
from subprocess import Popen, PIPE


def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="infile",
					  help="input file", metavar="INFILE")
	parser.add_option("-n", "--numlines", dest="numlines",
					  help="number of lines to generate (approx)", metavar="NUMLINES")
	(options, args) = parser.parse_args()
	return options

options = parse_options()
numlines = float(options.numlines)


p = Popen(['wc', '-l', options.infile], stdin=PIPE, stdout=PIPE, stderr=PIPE)
output, err = p.communicate()
line_count = output.split()[0]
line_count = int(line_count)
cutoff = numlines / line_count

f1 = options.infile
if (f1[-2:] == 'gz'):
	infile = gzip.open(f1, 'rt')
else:
	infile = open(f1, 'r')

for line in infile:
	line = line.rstrip()
	rand = random.random()
	if (rand < cutoff):
		print(line)

