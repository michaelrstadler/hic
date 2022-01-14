import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import gzip


############################################################################
# Utility functions
############################################################################

def load_viewer_file(infilename, norm=True, size=400):
        """Load viewer file of a Hi-C matrix, replace 0s and NaNs with dummy values, 
        trim outliers, take log.
        
        Args:
            infilename: string
                HiC matrix (panel) file
        
            size: int
                Expected size (rows and columns) of matrix
                
        Returns:
            x: numpy ndarray
                2D matrix of processed Hi-C data
        """
        # Load data.
        x = np.genfromtxt(infilename)
        # Resolve 0s and NaNs by setting to dummy value 0.1.
        dummy_val = 0.1
        x[x == 0] = dummy_val
        x[np.isnan(x)] = dummy_val
        # Remove large outlier values using a percentile cutoff.
        upper_trim = 99.9
        max_val = np.percentile(x, upper_trim)
        x[x > max_val] = max_val
        # Take the log.
        x = np.log(x)
        # Scale to a range of 0-1000.
        if (norm):
            x = (x - np.min(x)) / (np.max(x) - np.min(x)) * 1000
        return(x)

############################################################################
def load_track_data(trackfile_path):
    """Load genomic track data.

    Args:
        trackfile_path: str
            Path to genomic track file, format chr [tab] bin [tab] value

    Returns:
        track_data: dict of np arrays
            Keys are chromosomes, arrays contain values in bins matching
            size of input data
    """
    track_data = {}
    max_bin = {}
    with gzip.open(trackfile_path, 'rt') as infile:
        for line in infile:
            items = line.split()
            (chr_, bin_, val) = items
            chr_ = re.sub('chr', '', chr_)
            if (chr_ not in track_data):
                track_data[chr_] = np.zeros(int(1e8 / 500))
            bin_ = int(bin_)
            if (chr_ not in max_bin):
                max_bin[chr_] = 0
            if (bin_ > max_bin[chr_]):
                max_bin[chr_] = bin_
            if (bin_ < len(track_data[chr_])):
                track_data[chr_][bin_] = float(val)

        for chr_ in max_bin:
            track_data[chr_] = track_data[chr_][0:(max_bin[chr_]+1)]
    return track_data

############################################################################
def write_track_data(track_data, outfile):
    """Write genomic track data files from data.

    Args:
        track_data: dict of np arrays
            Keys are chromosomes, arrays contain values in bins
        outfile: str
            File to write to,format chr [tab] bin [tab] value.
            Only non-zero values are written.       
    """
    out = open(outfile, 'w')
    for chr_ in track_data:
        for i in range(0, len(track_data[chr_])):
            val = track_data[chr_][i]
            if (val > 0):
                out.write(chr_ + '\t' + str(i) + '\t' + str(val) + '\n')
    check_call(['gzip', outfile])

############################################################################
def panels_old_to_new_format(old_folder, new_folder):
    """Convert HiC panels from older format (with row and col names) to newer format.
    
    Args:
        old_folder: string
            Path to folder containing panel files in the old format
        new_folder: string
            Path of folder to write new files to
            
    Returns: None
    """
    old_files = os.listdir(old_folder)
    for old_file in old_files:
        if (re.search('.DS_Store', old_file)):
            continue
        x = np.genfromtxt(os.path.join(old_folder, old_file), skip_header=1)[:,1:]
        np.savetxt(os.path.join(new_folder, old_file),x,newline='\n', fmt='%.6e')

############################################################################
def load_panels(data_folder, binsize=500):
    """Load hi-c panels (viewer) data into memory.
    
    Args: 
        data_folder: folder
            Folder containing hi-c panels (viewer files)
        binsize: int
            Binsize (in base pairs) of the panels to load.
            
    Returns:
        panels: dict of dict of ndarray
            Keys are chromosome, then bin1 (left-most), value is
            2d numpy ndarray with hi-c data.
    """
    
    panels = {}
    for filename in os.listdir(data_folder):
        if (re.search(str(binsize) + '.txt', filename)):
            chr_, loc1, loc2, _ = filename.split('_')
            x = load_viewer_file(os.path.join(data_folder, filename))
            if (chr_ not in panels):
                panels[chr_] = {}
            panels[chr_][int(loc1)] = x
    return panels

############################################################################
def load_gff3(gff_filename):
    """Load gff3 file, storing chrom, position, score in pandas df
    
    Args:
        gff_filename: string
            Path to gff3 file
    
    Returns:
        gff: pandas DF
            Each row is gff entry, col 0: chromosome, 1: mean position,
            2: score
    """
    chr_col = 0
    posL_col = 3
    posR_col = 4
    score_col = 5
    gff_file = open(gff_filename, 'r')
    chrs = []
    pos = []
    scores = []
    
    for l in gff_file:
        if (not l[0] == '#'):
            items = l.split('\t')
            chr_ = items[chr_col]
            if (chr_[0:3] == 'chr'):
                chr_ = chr_[3:]
            if (chr_ in ['2L', '2R', '3L', '3R', 'X']):
                posL = int(items[posL_col])
                posR = int(items[posR_col])
                pos.append(np.mean([posL, posR]))
                chrs.append(chr_)
                scores.append(float(items[score_col]))
    gff = pd.DataFrame([chrs, pos, scores])
    return gff.T

############################################################################
def mesh_like(arr, n):
    """Make mesh grid for last n dimensions of an array
    
    Makes a meshgrid with the same shape as the last n dimensions of input
    array-like object.
    
    Args:
        arr: array-like
            Array-like object that has a shape parameter
        n: int
            Number of dimensions, from the right, for which to make meshgrid.
    
    Returns:
        meshes: list of ndarrays
            Each element of list corresponds to ordered dimension of input,
            ndarrays are corresponding meshgrids of same shape as arr.
    """
    if (n > arr.ndim):
        raise ValueError('n is larger than the dimension of the array')
    # Make vectors of linear ranges for each dimension.
    vectors = []
    for i in reversed(range(1, n+1)):
        a = np.arange(0, arr.shape[-i])
        vectors.append(list(a))
    # Make meshgrids from vectors.
    meshes = np.meshgrid(*vectors, sparse=False, indexing='ij')
    return meshes