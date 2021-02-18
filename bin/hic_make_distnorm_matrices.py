'''
From a folder of Hi-C viewer files, creates a distance normalization matrix
for each resolution. This is used to perform distance correction within
the Hi-C viewer.
'''
############################################################################
from optparse import OptionParser
import sys
import re
import numpy as np
import os
import sys
from importlib import reload
sys.path.append('/Users/michaelstadler/Bioinformatics/Projects/insulators/bin')
import hic_jupyter as hc

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--folder", dest="data_folder",
					  help="Folder containing viewer files", metavar="DATA_FOLDER")
	parser.add_option("-s", "--size", dest="size",
					  help="Size, in bins, of input data (each dimension)", metavar="SIZE",
					  default=400)

	(options, args) = parser.parse_args()
	return options

def get_bin_sizes(folder):
    """Find all bin sizes in folder."""
    bin_sizes = {}
    for f in os.listdir(data_folder):
        f = f.rstrip()
        items = f.split('_')
        if (len(items) == 4):
            binsize = int(re.sub('.txt.gz', '', items[3]))
            bin_sizes[binsize] = 1
    return bin_sizes.keys()


def make_distance_normalize_matrix(data_folder, size, bin_size):
    """Create distnorm matrix for a bin size."""
    binfile_match = '_' + str(bin_size) + '.txt'
    sums_ = np.zeros(size)
    counts = np.zeros(size)
    # Create a matrix of the distance to the diagonal.
    mesh = hc.mesh_like(np.zeros((size,size)), 2)
    diag_dist = abs(mesh[0] - mesh[1])

    for filename in os.listdir(data_folder):
        if binfile_match in filename:
            x = hc.load_viewer_file(os.path.join(data_folder, filename), norm=False)
            # For each distance off diagonal.
            for d in range(0, size):
            	# Get all values of this file at distance d from diagonal.
                vals_d = x[diag_dist == d]
                # Add the sum of these values to the running sum.
                sums_[d] = sums_[d] + np.sum(vals_d)
                # Add count of values to running count.
                counts[d] = counts[d] + len(vals_d)
    
    print('Done reading')
    
    # Get mean values.
    means = sums_ / counts
    
    # Build normalization matrix from means and diagonal distance matrix.
    norm_mat = np.zeros((size, size))
    for n in range(0, size):
        norm_mat[diag_dist == n] = means[n]
    # Normalize to 0-1000 (to match normalization of viewer files).
    norm_mat = (norm_mat - np.min(norm_mat)) / (np.max(norm_mat) - np.min(norm_mat)) * 1000
    # Write file.
    outfilepath = os.path.join(data_folder, 'distnorm_matrix_' + str(bin_size) + '.txt.gz')   
    np.savetxt(outfilepath, norm_mat, newline='\n', fmt='%.6e')

options = parse_options()
data_folder = options.data_folder
size = int(options.size)
bin_sizes = get_bin_sizes(data_folder)

for bin_size in bin_sizes:
    print(str(bin_size) + ' bp')
    make_distance_normalize_matrix(data_folder=data_folder, size=size, bin_size=bin_size)

