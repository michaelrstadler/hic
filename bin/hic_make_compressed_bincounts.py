#!/usr/bin/env python

"""
This script takes a paired alignment file and, assigns each end to a bin (some chunk of 
a chr_omosome defined by supplied bin size), and prints out the bin-bin counts for only
contacts within some width of the diagonal (distance between the bins). 

Prints a unique format. The file starts with a series of lines that start with a #.
These are the total bin counts, used to do normalization, if desired, subsequently. After 
that, the format is compressed:
chr_omosome    bin1    bin2    count

Only non-zero bin-bin linkages are printed.
"""

__author__      = "Michael Stadler"
__version__ = '1.0.0'

from optparse import OptionParser
import sys
import re
import gzip
import numpy as np
from scipy.sparse import dok_matrix

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--files", dest="filenames",
                      help="paired alignment files, comma separated", metavar="FILE")
    
    parser.add_option("-b", "--bin_size",
                       dest="bin_size", default=1000000,
                      help="bin size")
    
    parser.add_option("-w", "--width",
                       dest="width", default=1000,
                      help="width in bins from diagonal")

    
    parser.add_option("-s", "--file_stem",
                        dest="file_stem", default='none',
                      help="output file stem. Adds diag_bin_counts bin size and chr_")

    parser.add_option("-v", "--valid_chrs",
                        dest="valid_chrs", default=None,
                      help="Valid chr_omosome names, delimited by |")

    (options, args) = parser.parse_args()
    return options


def Add_read (chr_, bin1, bin2, num_bins):
    """Add a new read to bin_bin_counts."""
    if (chr_ not in bin_bin_counts):
        bin_bin_counts[chr_] = dok_matrix((num_bins, (2 * width) + 1))
    # Assign bin2 as a relative position to bin1.
    bin2_rel = (bin2 - bin1) + width
    bin_bin_counts[chr_][bin1, bin2_rel] = bin_bin_counts[chr_][bin1, bin2_rel] + 1
    

def add_to_totals(chr_, bin, num_bins):
    """Add count to total bin counts."""
    if (chr_ in bin_totals):
        bin_totals[chr_][bin] = bin_totals[chr_][bin] + 1
    else:
        bin_totals[chr_] = np.zeros(num_bins)

def update_max_bin(chr_, bin1, bin2):
    """Check if new read contains maximum observed bin for chr_omosome, update if so."""
    if (chr_ in max_bin):
        return(max(bin1, bin2, max_bin[chr_]))
    else:
        return(max(bin1, bin2))

# Main:
    
options = parse_options()
valid_chrs = options.valid_chrs
if valid_chrs is None:
    valid_chrs = []
else:
    valid_chrs = valid_chrs.split('|')

bin_size = int(options.bin_size)
filenames = options.filenames
files = filenames.split(',')
width = int(options.width)
num_bins = int(1e5) #dummy variable for array creation...chr_omosome entrants must be <50 million bp
bin_bin_counts = {}
bin_totals = {}
max_bin = {}

line_count = 0 # Keep track of lines processed.

# Read mapping file into memory.
print('Reading files:')
for f in files:
    print(f)
    if (f[-2:] == 'gz'):
        infile = gzip.open(f, 'rt')
    else:
        infile = open(f, 'r')

    for line in infile:
        line_count = line_count + 1
        if (line_count % 10000000 == 0):
            print('. ' + str(line_count / 10000000))
        line = line.rstrip()
        items = line.split()
        (chr1, Lmost1, chr2, Lmost2) = items[2], int(items[3]), items[5], int(items[6])
        if ((chr1 == chr2)):
            if ((len(valid_chrs) == 0) or (chr1 in valid_chrs)):
                bin1 = int(Lmost1 / bin_size)
                bin2 = int(Lmost2 / bin_size)
                add_to_totals(chr1, bin1, num_bins)
                if (bin1 != bin2):
                    add_to_totals(chr1, bin2, num_bins) 
                max_bin[chr1] = update_max_bin(chr1, bin1, bin2)
                if (abs(bin1 - bin2) <= width):
                    Add_read(chr1, bin1, bin2, num_bins) # Avoid double counting on diagonal. This will be triggered unless chr_omosome and bin are the same
                    if (bin1 != bin2):
                        Add_read(chr1, bin2, bin1, num_bins)
    infile.close()
print('done reading\n')

# Set up output file stem
file_stem = ''
if (options.file_stem == 'none'):
    file_stem = re.sub('.txt', '', files[0])
else:
    file_stem = options.file_stem

# Write output data to separate files for each chr_omosome.
print('Writing files:')
for chr_ in bin_totals.keys():
    print (chr_)
    with open(file_stem + '_CompressedBinCounts_' + str(bin_size) + 'bp_' + str(chr_) + '.txt','w') as outfile:
        # Print total bin counts.
        for bin in range(0, max_bin[chr_] + 1):
            if bin_totals[chr_][bin] > 0:
                outfile.write('#' + chr_ + '\t' + str(bin) + '\t')
                outfile.write(str(bin_totals[chr_][bin]) + '\n')

        # Print bin-bin counts.
        for bin1 in range(0, max_bin[chr_] + 1):
            for bin2 in range(0, (2 * width) + 1):
                bin2_orig = bin1 + (bin2 - width)
                if ((bin2_orig >= 0) and (bin2_orig <= max_bin[chr_])):
                    if (bin_bin_counts[chr_][bin1, bin2] > 0):
                        outfile.write(chr_ + '\t' + str(bin1) + '\t' + str(bin2_orig) + '\t' + str(bin_bin_counts[chr_][bin1, bin2]) + '\n')
    del(bin_bin_counts[chr_])
print("Done.")
