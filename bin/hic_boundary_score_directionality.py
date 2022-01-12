'''
Compute a boundary score from Hi-C matrices. 
'''
from optparse import OptionParser
import gzip
import numpy as np
import pandas as pd

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Reduced bin file", metavar="FILE")
    parser.add_option("-b", "--binsize", dest="bin_size", default=500,
                      help="Bin size in bp", metavar="BINSIZE")
    parser.add_option("-w", "--window_size", dest="window_size",
                    default=7, help="Window size for mean-differencing",
                    metavar="WINDOWSIZE")
    parser.add_option("-s", "--step_size", dest="step_size",
                    default=3, help="Step size for differencing",
                    metavar="STEPSIZE")
    parser.add_option("-c", "--coefficient", dest="coefficient",
                    default=2, help="Step size for differencing",
                    metavar="COEFFICIENT")
    parser.add_option("-n", action="store_true", dest="normalize",
                    help="Use normalized data", metavar="NORMALIZE")
    parser.add_option("-o", "--outstem", dest="outstem",
                    help="outfile stem", metavar="OUTSTEM")

    (options, args) = parser.parse_args()
    return options

options = parse_options()
outfilepath = options.outstem + '_insulationscore.txt'
window_size = int(options.window_size)
step_size = int(options.step_size)
shift_size = -1 * int(step_size / 2)
COEFFICIENT = float(options.coefficient)
NORMALIZE = options.normalize

def read_compressed_binfile(filepath):
    left = 0
    right = 1
    left_right_counts = np.zeros((int(1e5), 2))
    totals = np.zeros(int(1e5))
    max_bin = 0

    with gzip.open(filepath, 'rt') as file:
        for line in file:
            line = line.rstrip()
            # Handle total rows.
            if line[0] == '#':
                chr_, bin_, count = line[1:].split('\t')
                bin_ = int(bin_)
                totals[bin_] = float(count)
                if bin_ > max_bin:
                    max_bin = bin_

            else:
                chr_, bin1, bin2, count = line.split('\t')
                bin1, bin2 = [int(x) for x in [bin1, bin2]]
                count = float(count)
                if NORMALIZE:
                    count = (count * 10_000) / (totals[bin1] * totals[bin2])
                if bin1 > bin2:
                    left_right_counts[bin1, left] += count
                elif bin1 < bin2:
                    left_right_counts[bin1, right] += count
            
        directionality = np.log((left_right_counts[:,right] / (left_right_counts[:,left] + 0.001)) + 0.001)
        scale = (np.apply_along_axis(np.min, 1, left_right_counts) / np.mean(left_right_counts[:,right])) ** COEFFICIENT
        directionality_scaled = directionality * scale
        return directionality_scaled[:(max_bin+1)], chr_
    
directionality, chromosome = read_compressed_binfile(options.filename)
directionality = pd.Series(directionality)

insulation = directionality.rolling(window_size, center=True).mean().diff(step_size).shift(shift_size).to_numpy()
with open(outfilepath, 'w') as outfile:
    for n in range(len(insulation)):
        outfile.write(chromosome + '\t' + str(n) + '\t' + str(insulation[n]) + '\n')




