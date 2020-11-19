#!/usr/bin/env python
############################################################################

"""
Make viewer-ready panels of Hi-C matrixes from a compressed bin count file.


"""

__author__      = "Michael Stadler"
__version__ = '1.0.1'

from optparse import OptionParser
import sys
import re
import gzip
import os
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-w", "--windowsize", dest="windowsize", default=400,
					  help="Panel width, in bins", metavar="WINDOWSIZE")
	parser.add_option("-s", "--stepsize", dest="stepsize", default=200,
					  help="Step size, in bins", metavar="STEPSIZE")
	parser.add_option("-o", "--outfile_stem", dest="outfile_stem",
					  help="Path and chromosome for outfiles", metavar="OUTSTEM")
	(options, args) = parser.parse_args()
	return options

def get_max_bin_chr(file):
	"""Retrieve the chromosome name and maximum bin from input file."""
	with open(file, 'rb') as f:
	    f.seek(-2, os.SEEK_END)
	    while f.read(1) != b'\n':
	        f.seek(-2, os.SEEK_CUR) 
	    last_line = (f.readline().decode())
	last_line_items = last_line.split()
	return((last_line_items[0], int(last_line_items[1])))            	

def rebin_2d(arr, new_shape):
	"""Bin a 2D numpy array to fit a new shape, with each entry
	in the binned array equal to the sum of the constituent
	entries in the original."""
	shape = (new_shape[0], arr.shape[0] // new_shape[0],
		new_shape[1], arr.shape[1] // new_shape[1])
	return arr.reshape(shape).sum(-1).sum(1)

def rebin_1d(arr, new_shape):
    """Bin a 1D numpy array, with each entry in the binned array equal to 
    the sum of the constituent entries in the original"""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],)
    return arr.reshape(shape).sum(-1)

# Main:

options = parse_options()
input_binsize = 500
arb_constant = 10000000 #gets multiplied by the normalized coverage to avoid small number issues
bin_combinations = [1,2,4,8]
stepsize = int(options.stepsize)
windowsize = int(options.windowsize)
outfile_stem = options.outfile_stem
f = options.filename
(chromosome, max_bin) = get_max_bin_chr(f)
total_counts = np.zeros(max_bin + 1)
bin_bin_counts = np.zeros((max_bin+1, max_bin+1))

# Handle gzipped files.
if (f[-2:] == 'gz'):
	file1 = gzip.open(f, 'rt')
else:
	file1 = open(f, 'r')

# Read coverage file into memory.
for line in file1:
	line = line.rstrip()
	# Add #-marked lines to total bin counts.
	if (line[0] == '#'):
		line = line[1:]
		(chr, bin_, count) = line.split('\t')
		bin_ = int(bin_)
		# Very rarely you can have a bin show up in the # total counts lines that doesn't show up in the bin-bin linkages.
		# This occurs when the only contacts with this bin occur at a greater distance than the width. Because the size
		# of total counts is determined from the largest bin in bin-bin counts, this creates an error.
		if (bin_ < len(total_counts)):
			total_counts[int(bin_)] = int(float(count)) #int of '0.0' throws error...have to go through float for some reason
	# Add non-# lines to bin-bin counts.
	else:
		(chr, bin1, bin2, count) = line.split('\t')
		bin1 = int(bin1)
		bin2 = int(bin2)
		bin_bin_counts[bin1][bin2] =  int(float(count)) #int of '0.0' throws error...have to go through float for some reason
file1.close()	
print('Done reading file.')

# Bin, normalize, and write files for different bin combinations.
for bins_combine in bin_combinations:
	# using consistent terminology where 'lnew' extension means bin numbers for this particular size (combination), whereas
	# the same variables without 'new' or with 'orig' are in the bin numbers of the input data.
	stepsize_new = stepsize * bins_combine
	windowsize_new = windowsize * bins_combine
	binsize_new = input_binsize * bins_combine

	# Strategy for making panels: For each panel, make two matrices, a matrix of the bin-bin count data and a normalization
	# matrix where each entry is the product of the total counts for each involved bin. To perform vanilla normalization, 
	# simply divide the data matrix by the normalization matrix (and multiply by a suitable constant).
	for panel_start in range(0, max_bin, stepsize_new):
		panel_end = panel_start + windowsize_new
		if (panel_end <= max_bin):
			# Slice panel from bin-bin data; will be a larger matrix for larger bin combinations.
			panel = bin_bin_counts[panel_start:panel_end, panel_start:panel_end]
			# Bin panel to match desired window size. Binning to fixed size with different sized inputs
			# produces larger bins for larger bin combinations.
			panel_binned = rebin_2d(panel, [windowsize, windowsize])
			# Slice and bin total counts vector to match panel.
			total_counts_vector = total_counts[panel_start:panel_end]
			total_counts_vector_binned = rebin_1d(total_counts_vector, [windowsize])
			# Take outer product of total counts vector to produce normalization matrix.
			# Entry i,j = total counts bin i X total counts bin j.
			norm_matrix = np.outer(total_counts_vector_binned, total_counts_vector_binned)#
			# Replace 0s in normalization matrix with dummy value of 1 to avoid divide by 0.
			norm_matrix[norm_matrix == 0] = 1
			# Normalize by dividing matrices.
			panel_binned_norm = panel_binned / norm_matrix * arb_constant
			# Write file for normalized matrix.
			pos1 = panel_start * input_binsize
			pos2 = pos1 + windowsize_new * input_binsize
			outfilename = outfile_stem + '_' + str(pos1) + '_' + str(pos2) + '_' + str(binsize_new) + '.txt.gz'
			np.savetxt(outfilename, panel_binned_norm, newline='\n', fmt='%.6e')

