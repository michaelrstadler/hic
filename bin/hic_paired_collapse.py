"""
 This script takes a paired alignment file and does one thing:
     1) Collapses identical reads
          
 The collapsing is pretty brute force: just store everything in a hash and toss out 
 duplicates. Only tosses reads where both ends map to the exact same place. Needs a 
 high memory machine, it goes without saying. If I ever have huge datasets, will need
 to re-write with some sort of sorting strategy.
"""

from optparse import OptionParser
import sys
import re
import gzip

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="paired alignment files, comma separated", metavar="FILES")
    parser.add_option("-o", "--out_stem", dest="out_stem",
                      help="Stem for output file", metavar="OUTSTEM")
    (options, args) = parser.parse_args()
    return options

options = parse_options()
total_reads = 0.0
dup_reads = 0.0
distance_limit = 600

alignments = {}

infiles = options.filename.split(',')
out_stem = options.out_stem
#infile_handle = open(infile_name, 'r')
outfile_name = out_stem + 'Collapsed.txt'
outfile_handle = open(outfile_name, 'w')

for f in infiles:
    if (f[-2:] == 'gz'):
        file_handle = gzip.open(f, 'rt')
    else:
        file_handle = open(f, 'r')

    for line in file_handle:
        total_reads += 1
        line = line.rstrip()
        data = line.split('\t')
        combo_string = '_'.join(data[1:7])
        if(not combo_string in alignments):
            outfile_handle.write(line + '\n')
            alignments[combo_string] = 1
        else:
            dup_reads += 1
    file_handle.close()
        
outfile_handle.close()

outstuff = str(dup_reads) + '/' + str(total_reads) + ' reads, or ' + str(dup_reads / total_reads * 100.0) + '%, are duplicates\n'

sys.stdout.write(outstuff)