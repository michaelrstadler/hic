'''
Combine multiple compressed sum files by summing bin counts.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="comma-separated compressed coverage files.", metavar="FILES")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="Output file", metavar="OUTFILE")

	(options, args) = parser.parse_args()
	return options
            
def add_line(l, data):
    """Add data from a line to existing data."""
    chr_ = l[0]
    if (chr_[0] == '#'):
        bin_, val = l[1:3]
        if (chr_ not in data):
            data[chr_] = {}
        if (bin_ not in data[chr_]):
        	data[chr_][bin_] = 0
        data[chr_][bin_] = data[chr_][bin_] + float(val)
    else:
        bin1, bin2, val = l[1:4]
        if (chr_ not in data):
            data[chr_] = {}
        if (bin1 not in data[chr_]):
            data[chr_][bin1] = {}
        if (bin2 not in data[chr_][bin1]):
        	data[chr_][bin1][bin2] = 0
        data[chr_][bin1][bin2] = data[chr_][bin1][bin2] + float(val)

def add_file(filename, data):
	"""Add a new file of data."""
	f = open(filename, 'r')
	for l in f:
	    l = l.strip()
	    l = l.split()
	    add_line(l, data)
	f.close()

# Main.

options = parse_options()

files = options.filenames.split(',')
data = {}
for f in files[:-1]:
	add_file(f, data)

f = open(files[-1], 'r')
out = open(options.outfile, 'w')
for l in f:
    l = l.strip()
    if (l[0] == '#'):
        (chr_, bin_, val) = l.split()
        valsum = float(val) + float(data[chr_][bin_])
        out.write(chr_ + '\t' + bin_ + '\t' + str(valsum) + '\n')

    else:
        (chr_, bin1, bin2, val) = l.split()
        valsum = float(val) + float(data[chr_][bin1][bin2])
        out.write(chr_ + '\t' + bin1 + '\t' + bin2 + '\t' + str(valsum) + '\n')

f.close()

