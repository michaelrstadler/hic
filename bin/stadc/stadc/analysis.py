import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import gzip

############################################################################
# Analysis functions
############################################################################

def aggregate_signal_gff_landmarks(panels, gff, rad=20, binsize=500):
    """Aggregate the Hi-C signal around genomic positions in a gff file
    
    Produces the sum of the input hi-c signal around the supplied genomic
    positions.
    
    Args:
        panels: dict of dict of ndarray
            Keys are chromosome, then bin1 (left-most), value is
            2d numpy ndarray with hi-c data.
        gff: pandas DF
            Each row is gff entry, col 0: chromosome, 1: mean position,
            2: score
        rad: int
            "Radius" (halfwidth in x and y) of the region to be summed
        binsize: int
            Binsize of data to use (500, 1000, 2000, or 4000)
        
    Returns:
        sum_: ndarray
            Square array of dimension 2*rad + 1 containing summed signal
    """
    stepsize = binsize * 200
    lw = (2*rad) + 1
    # Initialize sum array with lengt/width according to rad.
    sum_ = np.zeros((lw, lw))
            
    for i in range(0, len(gff.index)):
        chr_, pos, score = gff.iloc[i]
        if (chr_ in panels):
            panel_start = int(pos / stepsize) * stepsize
            # If start is off of this panel, move back one step.
            if ((pos - (rad * binsize)) < panel_start):
                panel_start = panel_start - stepsize
            if (panel_start >= 0 and panel_start in panels[chr_]):
                panel = panels[chr_][panel_start]
                bin_ = int(pos / binsize) - int(panel_start / binsize) + 1
                start = bin_ - rad
                stop = bin_ + rad + 1
                sub = panel[start:stop, start:stop]
                sub[np.isnan(sub)] = 0
                sum_ += sub

    return sum_  


############################################################################
def aggregate_signal_gff_landmarks_norm(panels, gff, rad=99, binsize=500):
    """aggregate_signal_gff_landmarks with normalization
    
    Wrapper for aggregate_signal_gff_landmarks. Generates a "randomized"
    normalization matrix by permuting chromosome labels. Randomized matrix
    accounts for distance decay (and potentially other features of data).
    
    Args:
        panels: dict of dict of ndarray
            Keys are chromosome, then bin1 (left-most), value is
            2d numpy ndarray with hi-c data.
        gff: pandas DF
            Each row is gff entry, col 0: chromosome, 1: mean position,
            2: score
        rad: int
            "Radius" (halfwidth in x and y) of the region to be summed
        binsize: int
            Binsize of data to use (500, 1000, 2000, or 4000)
        
    Returns:
        signal: ndarray
            Square array of dimension 2*rad + 1 containing summed signal
    
    """
    def swap_chrs_gff(gff_orig):
        """Make shuffled version of a gff data structure by 'shifting' each
        chromosome by one (not random--chromosome will never be itself)"""
        gff = gff_orig.copy()
        chrs = np.unique(gff[0])
        convert = {}
        convert[chrs[0]] = chrs[-1]
        for i in range(1,len(chrs)):
            convert[chrs[i]] = chrs[i-1]
        for i in range(0, len(gff.index)):
            gff.iloc[i,0] = convert[gff.iloc[i,0]]
        return gff

    signal_real = aggregate_signal_gff_landmarks(panels, gff, rad, binsize)
    signal_shuff = aggregate_signal_gff_landmarks(panels, swap_chrs_gff(gff), rad, binsize)
    # Subtract because inputs (normalized counts) are already log.
    return signal_real - signal_shuff

############################################################################
def downsample_match_compressed_bin(small_filename, big_filename, outfilename, width=3500):
    """Downsample a compressed_bin file to match a lower-coverage file
    
    Args:
        small_filename: string
            Compressed bin file with fewer counts
        big_filename: string
            Compressed bin file with more counts
        outfilename: string
            File to write down-sampled output to
        width: int
            Width, in bins, of input data
    
    Returns: None
    """
    
    def get_count(bin_bin_counts, bin1, bin2):
        """Get count for a single orientation of a bin-bin linkage"""
        if (bin1 in bin_bin_counts):
            if (bin2 in bin_bin_counts[bin1]):
                return bin_bin_counts[bin1][bin2]
        return 0
    
    # Get the total counts in the smaller file (note, bins will 
    # be double counted but since we are just looking for a comparison, it works)
    total_smallfile = 0
    f1 = open(small_filename, 'r')
    for line in f1:
        line = line.rstrip()
        if (line[0] == '#'):
            continue
        items = line.split()
        count = float(items[3])
        total_smallfile += count
    f1.close()

    # Read in data from larger file. bin_names is a list that stores the first three columns
    # of the file, counts is a list that stores the counts. They are index-matched, so that 
    # we can draw from counts and use the indexes to get the corresponding bins.
    bin_pair_names = []
    counts = np.array([])
    max_bin = 0
    f2 = open(big_filename, 'r')
    for line in f2:
        line = line.rstrip()
        if (line[0] == '#'):
            continue
        items = line.split()
        # Update max bin.
        max_ = max(int(items[1]), int(items[2]))
        if (max_ > max_bin):
            max_bin = max_
        bin_pair_name = '\t'.join(items[0:3])
        count = float(items[3])
        bin_pair_names.append(bin_pair_name)
        counts = np.append(counts,count)
        
    f2.close()
    
    # Draw reads from larger file.
    """
    There is a complication here. I'm drawing from the bin pairs independently, but
    they're not independent — bin 10-bin 14 is the same as bin 14-bin 10. To compensate
    for this, we divide the number of draws by 2. When we go to write the outputs, we
    will sum the draws from the two orientations of each bin pair.
    """
    total_reads = counts.sum()
    number_to_draw = int(total_smallfile/ 2)
    probs = counts / total_reads
    # Randomly sample indexes (0 to length of bin_pair_names), without replacement, according
    # to probabilities based on counts in the original file.
    draw = np.random.choice(np.arange(0,len(bin_pair_names)), number_to_draw,
              p=probs, replace=False)
    # Count up the occurences of each index.
    indexes, index_counts = np.unique(draw, return_counts=True)
    
    # Get chromosome name.
    chr_ = bin_pair_names[0].split('\t')[0]
    
    # From the sampled indexes and counts, build up bin_bin_counts 
    # and bin_counts to store paired linkage counts and total counts.
    bin_counts = np.zeros(max_bin + 1)
    bin_bin_counts = {}
    for i in range(0, len(indexes)):
        index = indexes[i]
        count = index_counts[i]
        bin_name = bin_pair_names[index]
        chr, bin1, bin2 = bin_name.split()
        bin1 = int(bin1)
        bin2 = int(bin2)
        
        if (bin1 not in bin_bin_counts):
            bin_bin_counts[bin1] = {}
        bin_bin_counts[bin1][bin2] = count

        bin_counts[bin1] += 1
        if (bin1 != bin2):
            bin_counts[bin2] += 1
    
    # Write output with total counts and bin linkage counts.
    outfile = open(outfilename, 'w')
    for i in range(0, len(bin_counts)):
        outfile.write('#' + chr_ + '\t' + str(i) + '\t' + str(bin_counts[i]) + '\n')
            
    for bin1 in range(0, max_bin):
        for bin2 in range(max(0, bin1 - width), min(max_bin, bin1 + width)):
            # Get counts for bin pair from each orientation.
            count1 = get_count(bin_bin_counts, bin1, bin2)
            count2 = get_count(bin_bin_counts, bin2, bin1)
            symmetric_count = count1 + count2
            if (symmetric_count > 0):
                outfile.write(chr_ + '\t' + str(bin1) + '\t' + str(bin2) + '\t' + str(symmetric_count) + '\n')
        
        
    outfile.close()

############################################################################
def insulation_cum_dist(agg_orig, buffer=2, min_dist=3, max_dist=100, binsize=500):
    """Calculate an insulation score as a function of distance from an aggregate 
    Hi-C matrix
    
    Description: Starts with a matrix of the aggregate Hi-C data in which the feature
    is in the center of the matrix. Insulation score compares the signal at a fixed 
    genomic distance either across the feature or in the opposite direction (agnostic
    to whether features present or not). Visually, one can think of the aggregate matrix
    in four quadrants, (I to IV clockwise from top left), with quadrants II and IV
    representing interactions across the feature and I and III being not across it. 
    To get the score at a given distance, we draw a line of slope -1 at a the desired 
    distance from the diagonal, and compare the mean signal of the pixels in the "across"
    region to the "not across" region.
    
    Args:
        agg_orig: numpy ndarray
            Aggregate Hi-C matrix
        buffer: int
            Number of bins on either side of the boundary (feature) to omit in calculation
        min_dist: int
            Starting distance from diagonal to compute
        max_dist: int
            Ending distance from diagonal to compute

    Returns:
        scores: list
            Insulation scores by distance
        distances: list
            Matching genomic distances
    """
    agg = agg_orig.copy()
    # Reset minimum value to 0 to avoid negative number problems.
    if (agg.min() < 0):
        agg = agg - agg.min()
    # Make mesh grid.
    y, x = mesh_like(agg, 2)
    # Make a grid of distances to the diagonal, with positive values to upper right.
    dist_diag = x-y
    midpoint = int(agg.shape[0] / 2)
    # Define insulated region as top right quadrant (when displayed as an image: y-axis is backwards).
    insulated_region = (x > (midpoint + buffer)) & (y < (midpoint - buffer))
    noninsulated_region = (x < (midpoint - buffer)) | (y > (midpoint + buffer))
    scores = []
    
    # Go through range of distances off diagonal, compute insulation at each distance.
    for dist in range(min_dist, max_dist):
        insulated_indexes = np.where((dist_diag == dist) & insulated_region)
        # Restrict non-insulated pixels to the regions where one of the points at 
        # distance dist is across the midpoint.
        xlim = midpoint - dist
        ylim = midpoint + dist
        non_insulated_indexes = np.where(noninsulated_region & (x >= xlim) & (y <= ylim) & (dist_diag == dist))
        insulated = agg[insulated_indexes]  
        non_insulated = agg[non_insulated_indexes]
        # Is division the right thing here?
        scores.append(non_insulated.mean() / insulated.mean())
    distances = np.arange(min_dist, max_dist) * binsize
    return scores, distances

############################################################################
def features_assign_chip_scores(feature_file, signal_files, outfolder, 
            outputname, names, norm=True, width=2, feature_file_cols=[0,1], 
            signal_file_cols=[0,1,3], max_chr_size = 5e7, binsize=500):
    """Assign relative scores to genomic features from ChIP-like input data
    
    The concept is that boundaries are sites of combinatoric insulator protein
    binding, but that different different boundaries have distinct combinations
    of protein binding, with individual proteins more "dominant" at some
    boundaries than others. To this end, input ChIP signals are first normalized,
    then each boundary (or other feature) is assigned a score for each protein
    as a fraction of the total ChIP signal at that boundary. So, protein X
    will have a high score at a boundary where its signal is high but there is 
    little signal for the other proteins in the set. The output is a gff-like 
    file for each protein with the fractional score in the score column. This 
    file can be used as input to other scripts and functions.
    
    Args:
        feature_file:
            File containing genomic features (boundaries)
        signal_files: list-like
            Files containing ChIP-like data (genome-wide, continuous)
        outfolder: string
            Path to folder to write outfile
        outputname: string
            Name for output files to which names entry is appended.
        names: list-like, iterable
            Names corresponding to the signal files to be used in output files
        norm: bool
            If true, reports ChIP signal as a fraction of total ChIP signal. 
            Otherwise, just reports the ChIP signal itself.
        width: int
            Width, in bins, of the region used to calculate the signal at each 
            feature
        feature_file_cols: list
            Columns in the feature file corresponding to [chromosome, position]
        signal_file_cols: list
            Columns in the signal files corresponding to [chromosome, position
            score/value]
        max_chr_size: int
            Max size of chromosomes, in base pairs (used for memory allocation)
        binsize: int
            Size, in base pairs, of bins in signal files
        
        
    """
    # Read in data from signal files.
    max_chr_bin = int(max_chr_size / binsize)
    # Store signal in dict of ndarray, keys are chromosomes, dim1 is signal ID (
    # one per signal file) dim2 is genomic position.
    signal = {}
    for i in range(0, len(signal_files)):
        signal_file = signal_files[i]
        with open(signal_file, 'r') as f:
            for line in f:
                chr_, loc, val = np.array(line.rstrip().split())[signal_file_cols]
                val = float(val)
                chr_ = re.sub('chr', '', chr_)
                if (chr_ not in signal):
                    signal[chr_] = np.zeros((max_chr_bin, len(signal_files)))
                bin_ = int(int(loc) / binsize)
                signal[chr_][bin_, i] = val
    
    # Normalize signals.
    for chr_ in signal.keys():
        sums = np.sum(signal[chr_], axis=0)
        sums_mean = np.mean(sums)

        signal[chr_] = signal[chr_] / sums * sums_mean
    
    # Read features into dict of list, keys are chromosomes, list entries are
    # genomic positions.
    features = {}
    with open(feature_file, 'r') as f:
        for line in f:
            if (line[0] == '#'):
                continue
            items = np.array(line.rstrip().split())
            chr_, loc = items[feature_file_cols]
            chr_ = re.sub('chr', '', chr_)
            bin_ = int(int(loc) / binsize)
            if (chr_ not in features):
                features[chr_] = []
            features[chr_].append(bin_)
    
    # Write out file for each input signal.
    for i in range(0, len(signal_files)):
        outfile = os.path.join(outfolder, outputname + names[i] + '.txt')
        with open(outfile, 'w') as out:
            for chr_ in features:
                for bin_ in features[chr_]:
                    # Take the mean value for each signal in the window defined by width.
                    means = np.mean(signal[chr_][(bin_-width):(bin_+width+1)], axis=0)
                    # For each signal, value is either the mean signal or its fraction 
                    # of the total signal.
                    if (norm):
                        value = means[i] / np.sum(means)
                    else:
                        value = means[i]
                    pos = str(bin_ * binsize)
                    out.write(chr_ + '\t.\t.\t' + pos + '\t' + pos + '\t' + str(value) + '\t.\t.\t.\n')
                    
############################################################################
def plot_aggregate_features(panels, gff, rad=99, binsize=500):
    """Plot aggregate Hi-C signal and insulation score for a genomic feature set
    
    Plots the aggregate (sum/mean) Hi-C score around a set of features supplied
    in a gff file, and uses the score feature of the gff to split the features
    into thirds and plots the same for these splits. The top row is for the
    entire set, row 2-4 are for thirds with successively higher scores.
    
    Args:
        panels: dict of dict of ndarray
            Hi-C data panels. Keys are chromosome, then bin1 (left-most), value
            is 2d numpy ndarray with hi-c data.
        gff: pandas DF
            gff object produced by load_gff. Each row is gff entry, col 0: 
            chromosome, 1: mean position, 2: score
        rad: int
            The 'radius' of the window for the aggregate score. Length of
            matrix will be 2*rad + 1.
        binsize: int
            The bin size, in bp, of the data to use (500, 1000, 2000, 4000)
    
    Returns: None
    """
    t1, t2 = np.percentile(list(gff.iloc[:,2]), (33, 66))
    gff_1 = gff[gff.iloc[:,2] <= t1]
    gff_2 = gff[(gff.iloc[:,2] > t1) & (gff.iloc[:,2] < t2)]
    gff_3 = gff[gff.iloc[:,2] >= t2]
    signal_norm_all = aggregate_signal_gff_landmarks_norm(panels, gff, rad, binsize)
    signal_norm_1 = aggregate_signal_gff_landmarks_norm(panels, gff_1, rad, binsize)
    signal_norm_2 = aggregate_signal_gff_landmarks_norm(panels, gff_2, rad, binsize)
    signal_norm_3 = aggregate_signal_gff_landmarks_norm(panels, gff_3, rad, binsize)
    
    fig, ax = plt.subplots(4,2, figsize=(15,15))
    ax[0][0].imshow(signal_norm_all)
    ax[1][0].imshow(signal_norm_1)
    ax[2][0].imshow(signal_norm_2)
    ax[3][0].imshow(signal_norm_3)
    ax[0][1].plot(insulation_cum_dist(signal_norm_all)[0])
    ax[1][1].plot(insulation_cum_dist(signal_norm_1)[0])
    ax[2][1].plot(insulation_cum_dist(signal_norm_2)[0])
    ax[3][1].plot(insulation_cum_dist(signal_norm_3)[0])

############################################################################
def plot_aggregate_features_compare(panels1, panels2, gff, rad=99, binsize=500):
    """Plot aggregate Hi-C signal and insulation score for a genomic feature set
    in two Hi-C datasets
    
    Plots the aggregate (sum/mean) Hi-C score around a set of features supplied
    in a gff file as well as the insulation-distance plot for the aggregate data.
    Uses the score feature of the gff to split the features into thirds and plots 
    the same for these splits. The top row is for the entire set, row 2-4 are for 
    thirds with successively higher scores. Aggregate plots for panels1 is plotted 
    in col 0, panels2 in col 1, and insulation-distance plots are overlaid in 
    col 2. Data for panels1 is in blue, for panels2 in orange.
    
    Args:
        panels1, panels2: dict of dict of ndarray
            Hi-C data panels. Keys are chromosome, then bin1 (left-most), value
            is 2d numpy ndarray with hi-c data.
        gff: pandas DF
            gff object produced by load_gff. Each row is gff entry, col 0: 
            chromosome, 1: mean position, 2: score
        rad: int
            The 'radius' of the window for the aggregate score. Length of
            matrix will be 2*rad + 1.
        binsize: int
            The bin size, in bp, of the data to use (500, 1000, 2000, 4000)
    
    Returns: None
    """
    t1, t2 = np.percentile(list(gff.iloc[:,2]), (33,66))
    gff_1 = gff[gff.iloc[:,2] <= t1]
    gff_2 = gff[(gff.iloc[:,2] > t1) & (gff.iloc[:,2] < t2)]
    gff_3 = gff[gff.iloc[:,2] >= t2]

    sig1_norm_all = aggregate_signal_gff_landmarks_norm(panels1, gff, rad, binsize)
    sig1_norm_1 = aggregate_signal_gff_landmarks_norm(panels1, gff_1, rad, binsize)
    sig1_norm_2 = aggregate_signal_gff_landmarks_norm(panels1, gff_2, rad, binsize)
    sig1_norm_3 = aggregate_signal_gff_landmarks_norm(panels1, gff_3, rad, binsize)
    sig2_norm_all = aggregate_signal_gff_landmarks_norm(panels2, gff, rad, binsize)
    sig2_norm_1 = aggregate_signal_gff_landmarks_norm(panels2, gff_1, rad, binsize)
    sig2_norm_2 = aggregate_signal_gff_landmarks_norm(panels2, gff_2, rad, binsize)
    sig2_norm_3 = aggregate_signal_gff_landmarks_norm(panels2, gff_3, rad, binsize)
    
    Rmost = (rad + 1) * binsize / 1000
    Lmost = -1 * Rmost
    min_dist = 10
    x_vals = np.arange(min_dist, rad+1) * binsize / 1000
    fig, ax = plt.subplots(4,3, figsize=(15,15))
    img1 = ax[0][0].imshow(sig1_norm_all, cmap='inferno')
    img1.set_extent([Lmost, Rmost, Lmost, Rmost])
    img2 = ax[1][0].imshow(sig1_norm_1, cmap='inferno')
    img2.set_extent([Lmost, Rmost, Lmost, Rmost])
    img3 = ax[2][0].imshow(sig1_norm_2, cmap='inferno')
    img3.set_extent([Lmost, Rmost, Lmost, Rmost])
    img4 = ax[3][0].imshow(sig1_norm_3, cmap='inferno')
    img4.set_extent([Lmost, Rmost, Lmost, Rmost])
    img5 = ax[0][1].imshow(sig2_norm_all, cmap='inferno')
    img5.set_extent([Lmost, Rmost, Lmost, Rmost])
    img6 = ax[1][1].imshow(sig2_norm_1, cmap='inferno')
    img6.set_extent([Lmost, Rmost, Lmost, Rmost])
    img7 = ax[2][1].imshow(sig2_norm_2, cmap='inferno')
    img7.set_extent([Lmost, Rmost, Lmost, Rmost])
    img8 = ax[3][1].imshow(sig2_norm_3, cmap='inferno')
    img8.set_extent([Lmost, Rmost, Lmost, Rmost])
    ax[0][2].plot(x_vals, insulation_cum_dist(sig1_norm_all, min_dist=min_dist)[0])
    ax[0][2].plot(x_vals, insulation_cum_dist(sig2_norm_all, min_dist=min_dist)[0])
    ax[1][2].plot(x_vals, insulation_cum_dist(sig1_norm_1, min_dist=min_dist)[0])
    ax[1][2].plot(x_vals, insulation_cum_dist(sig2_norm_1, min_dist=min_dist)[0])
    ax[2][2].plot(x_vals, insulation_cum_dist(sig1_norm_2, min_dist=min_dist)[0])
    ax[2][2].plot(x_vals, insulation_cum_dist(sig2_norm_2, min_dist=min_dist)[0])
    ax[3][2].plot(x_vals, insulation_cum_dist(sig1_norm_3, min_dist=min_dist)[0])
    ax[3][2].plot(x_vals, insulation_cum_dist(sig2_norm_3, min_dist=min_dist)[0])