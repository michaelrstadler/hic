#!/usr/bin/env python

"""
This script takes a paired alignment file and, assigns each end to a bin (some chunk of 
a chr_omosome defined by supplied bin size). Reads are stored in a sparse matrix format
and saved.

Two types of files are saved:

    1. bin_totals.txt stores total counts for bins of all chromosomes, used for 
        normalization. Format is:
            chromosome    bin    count
    
    2. Sparse matrixes for each chromosome. Files are gzipped pickled scipy.sparse dok_matrix. 
        Filenames are in format: 2R_sparsemat.pkl.gz
"""

__author__      = "Michael Stadler"
__version__ = '1.0.0'

from optparse import OptionParser
import sys
import os
import re
import gzip
import numpy as np
import pickle
from time import time
from scipy.sparse import dok_matrix

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--files", dest="filenames",
                      help="paired alignment files, comma separated", metavar="FILE")
    
    parser.add_option("-b", "--bin_size",
                       dest="bin_size", default=1000000,
                      help="bin size")
    
    parser.add_option("-o", "--out_folder",
                        dest="out_folder",
                      help="Folder for output")

    parser.add_option("-v", "--valid_chrs",
                        dest="valid_chrs", default=None,
                      help="[optiona]Valid chr_omosome names, delimited by |")

    (options, args) = parser.parse_args()
    return options


def Add_read (chr_, bin1, bin2, num_bins):
    """Add a new read to bin_bin_counts."""
    if (chr_ not in bin_bin_counts):
        bin_bin_counts[chr_] = dok_matrix((num_bins, num_bins))
    # Assign bin2 as a relative position to bin1.
    bin_bin_counts[chr_][bin1, bin2] = bin_bin_counts[chr_][bin1, bin2] + 1
    

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
t1 = time()
options = parse_options()

out_folder = options.out_folder
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)

valid_chrs = options.valid_chrs
if valid_chrs is None:
    valid_chrs = []
else:
    valid_chrs = valid_chrs.split('|')

bin_size = int(options.bin_size)
filenames = options.filenames
files = filenames.split(',')
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
        if (line_count % 1000000 == 0):
            sys.stdout.write('.')
            sys.stdout.flush()

        line = line.rstrip()
        items = line.split()
        (chr1, Lmost1, chr2, Lmost2) = items[2], int(items[3]), items[5], int(items[6])
        if (chr1 == chr2):
            # If valid_chrs is empty, all are taken. Otherwise, must be in valid_chrs.
            if ((len(valid_chrs) == 0) or (chr1 in valid_chrs)):
                bin1 = int(Lmost1 / bin_size)
                bin2 = int(Lmost2 / bin_size)
                max_bin[chr1] = update_max_bin(chr1, bin1, bin2)
                # Add to bin totals.
                add_to_totals(chr1, bin1, num_bins)
                # Avoid double counting diagonal.
                if (bin1 != bin2):
                    add_to_totals(chr1, bin2, num_bins) 
                # Add to sparse matrices.
                Add_read(chr1, bin1, bin2, num_bins) 
                # Avoid double counting on diagonal.
                if (bin1 != bin2):
                    Add_read(chr1, bin2, bin1, num_bins)
    infile.close()
print('done reading\n')

# Write bin totals file.
bin_totals_file = os.path.join(out_folder, 'bin_totals.txt')
with open(bin_totals_file, 'w') as outfile:
    for chr_ in bin_totals:
        for bin in range(0, max_bin[chr_] + 1):
            if bin_totals[chr_][bin] > 0:
                outfile.write(chr_ + '\t' + str(bin) + '\t')
                outfile.write(str(bin_totals[chr_][bin]) + '\n')

# Write sparse matrix files for individual chromosomes.
for chr_ in bin_bin_counts:
    chr_file = os.path.join(out_folder, chr_ + '_sparsemat.pkl.gz')
    with gzip.open(chr_file, 'wb') as file:
        pickle.dump(bin_bin_counts[chr_], file, protocol=4)

t2 = time()
print('time: ' + str(t2 - t1))
