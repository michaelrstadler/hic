#!/usr/bin/env python
############################################################################

"""
Mask Hi-C panels based on a file of masked regions of the genome that have
low visibility. Masked positions are set to np.nan.

Inputs:
    input_folder: 
        Folder containing viewer panel files. Hi-C panels are masked,
        all other files are copied to output folder as is.
    output_folder: 
        Folder to write masked files to
    mask_file:
        File containing regions to mask. Line format is:
        chromosome [tab] bin number (500 bp) [tab] 1

"""
__author__      = "Michael Stadler"
__version__ = '1.0.0'

from optparse import OptionParser
import gzip
import re
import os
from shutil import copyfile
import numpy as np

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--infolder", dest="input_folder",
                      help="Folder containing viewer files", metavar="INFOLDER")
    parser.add_option("-o", "--outfolder", dest="output_folder",
                      help="Folder to write new files to", metavar="OUTFOLDER")
    parser.add_option("-m", "--maskfile", dest="mask_file",
                      help="File containing masked regions", metavar="MASKFILE")
    (options, args) = parser.parse_args()
    return options

def load_track_data(trackfile_path):
    """Load genomic track data."""
    track_binsize = 500
    track_data = {}
    with gzip.open(trackfile_path, 'rt') as infile:
        for line in infile:
            items = line.split()
            (chr_, bin_, val) = items
            chr_ = re.sub('chr', '', chr_)
            if (chr_ not in track_data):
                track_data[chr_] = np.zeros(int(1e8 / 500))
            bin_ = int(bin_)
            if (bin_ < len(track_data[chr_])):
                track_data[chr_][bin_] = float(val)
    return track_data

def rebin_1d(arr, new_shape):
    """Bin a 1D numpy array, with each entry in the binned array equal to 
    the mean of the constituent entries in the original"""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],)
    return arr.reshape(shape).mean(-1)

options = parse_options()
mask_file = options.mask_file
input_folder = options.input_folder
output_folder = options.output_folder

# Bin the 500 bp mask to form other resolutions, using a threshold of half of the bin
# being masked to mask the larger bin.
mask_500 = load_track_data(mask_file)
mask = {}
mask[500] = mask_500
binsizes = [500, 1000, 2000, 4000]
for i in range(1,4):
    new_binsize = binsizes[i]
    mask[new_binsize] = {}
    bin_factor = 2 ** i
    for chr_ in mask_500.keys():
        new_size = int(len(mask_500[chr_]) / bin_factor)
        binned = rebin_1d(mask_500[chr_], [new_size])
        mask[new_binsize][chr_] = np.where(binned >= 0.5, 1, 0)

# Make output folder if doesn't exist. 
try: 
    os.mkdir(output_folder) 
except: 
    pass

# Mask all viewer panel files, copy other files.  
files = os.listdir(input_folder)
for f in files:
    items = f.split('_')
    # If format matches viewer panel file.
    if (len(items) == 4):
        chr_, pos1, pos2, _ = items
        binsize = int(_.split('.')[0])
        bin1 = int(int(pos1) / binsize)
        bin2 = int(int(pos2) / binsize)
        x = np.genfromtxt(os.path.join(input_folder, f))
        # Get indexes of masked bins within this panel.
        masked = np.where(mask[binsize][chr_][bin1:bin2] > 0)
        # Set masked columns and rows to nan.
        x[masked, :] = np.nan
        x[:, masked] = np.nan
        outfilepath = os.path.join(output_folder, f)
        np.savetxt(outfilepath, x, newline='\n', fmt='%.6e')
    else:
        old = os.path.join(input_folder, f)
        new = os.path.join(output_folder, f)
        copyfile(old, new)
        