#!/usr/bin/env python
############################################################################

"""
Make viewer-ready panels of Hi-C matrixes from sparse matrix files. Viewer
panels are text forms of numpy ndarrays with chromosome and position 
indicated in filenames. Files are normalized via the vanilla normalization 
style.
"""

__author__      = "Michael Stadler"
__version__ = '1.1.0'

from optparse import OptionParser
import pickle
import sys
import re
import gzip
import os
import scipy
import numpy as np

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--folder", dest="infolder",
                      help="Folder with sparse matrix files", metavar="FOLDER")
    parser.add_option("-o", "--outfolder", dest="outfolder",
                      help="Path for output files", metavar="OUTFOLDER")
    parser.add_option("-w", "--windowsize", dest="windowsize", default=400,
                      help="Panel width, in bins", metavar="WINDOWSIZE")
    parser.add_option("-s", "--stepsize", dest="stepsize", default=200,
                      help="Step size, in bins", metavar="STEPSIZE")
    
    (options, args) = parser.parse_args()
    return options              

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

def read_bin_totals(bin_totals_file):
    """Read bin total from file, return in dict."""
    total_counts = {}
    max_bin = {}
    with open(bin_totals_file, 'r') as file:
        for line in file:
            chr_, bin_, count = line.rstrip().split('\t')
            bin_ = int(bin_)
            count = float(count)
            if chr_ in total_counts:
                total_counts[chr_][bin_] += count
                if bin_ > max_bin[chr_]:
                    max_bin[chr_] = bin_
            else:
                max_bin[chr_] = bin_
                total_counts[chr_] = np.zeros(int(1e5))
                total_counts[chr_][bin_] = count
    return total_counts, max_bin

def pad(panel, windowsize):
    """Pad an array with zeros if it goes over the edge of the matrix."""
    padded = np.pad(panel, 
        ([0, windowsize - panel.shape[0]],[0, windowsize - panel.shape[1]])
        )
    return padded

####
####

if __name__ == "__main__":

    options = parse_options()
    input_binsize = 500
    arb_constant = 10000000 #gets multiplied by the normalized coverage to avoid small number issues
    bin_combinations = [1,2,4,8, 16, 32, 64, 128]
    stepsize = int(options.stepsize)
    windowsize = int(options.windowsize)
    outfolder = options.outfolder
    infolder = options.infolder

    # Read total counts file.
    bin_totals_file = os.path.join(infolder, 'bin_totals.txt')
    total_counts, max_bin = read_bin_totals(bin_totals_file)

    # Normalize and write viewer panel files.
    for file in os.listdir(infolder):
        if file[-17:] == '_sparsemat.pkl.gz':
            print(file)
            chr_ = file[:-17]
            filepath = os.path.join(infolder, file)
            with gzip.open(filepath, 'rb') as pf:
                mat = pickle.load(pf)

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
                for panel_start in range(0, max_bin[chr_], stepsize_new):
                    panel_end = panel_start + windowsize_new
                    # CHANGE CHANGE CHANGE!!!!!!!!!!!!!!!!
                    if (panel_start <= max_bin[chr_]):
                        # Slice panel from bin-bin data; will be a larger matrix for larger bin combinations.
                        panel = mat[panel_start:panel_end, panel_start:panel_end]
                        # Convert sparse matrix to numpy array.
                        
                        panel = panel.toarray()
                        # Bin panel to match desired window size. Binning to fixed size with different sized inputs
                        # produces larger bins for larger bin combinations.
                        if panel.shape != (windowsize_new, windowsize_new):
                            panel = pad(panel, windowsize_new)
                        panel_binned = rebin_2d(panel, [windowsize, windowsize])
                        # Slice and bin total counts vector to match panel.
                        total_counts_vector = total_counts[chr_][panel_start:panel_end]
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
                        outfilename = os.path.join(outfolder , chr_ + '_' + str(pos1) + '_' + str(pos2) + '_' + str(binsize_new) + '.txt.gz')
                        np.savetxt(outfilename, panel_binned_norm, newline='\n', fmt='%.6e')
