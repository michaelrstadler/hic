#!/usr/bin/env python

"""sort_GArepeats.py: Search paired fastq files for pairs that have GA
repeats on one read, write two fastq files containing corresponding
reads.


TO DO:
-
"""

__author__      = "Michael Stadler"
__copyright__   = "Copyright 2022, California, USA"
__version__ = "1.0.0"

import argparse
import re
import os
import gzip
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument("-f", "--file1", type=str,  required=True,
                help="First fastq file.")
    parser.add_argument("-s", "--file2", type=str,  required=True,
                help="Second fastq file.")
    parser.add_argument("-o", "--outdir", type=str,  required=True,
                help="Directory to which to write output fastq files.")
    parser.add_argument("-p", "--pattern", type=str,  required=False,
                default="(AAGAG){5}", 
                help="Directory to which to write output fastq files.")
    
    args = parser.parse_args()
    return args

args = parse_args()

# Check inputs.
for f in [args.file1, args.file2]:
    if f[-8:] != 'fastq.gz':
        raise ValueError('Input files must be .fastq.gz')   

    if not os.path.exists(f):
        raise ValueError('File ' + f + ' does not exist.')
    
if not os.path.isdir(args.outdir):
    raise ValueError('Outdir is not a valid directory.')


# Setup outfiles.
stem1 = re.sub('.gz', '', args.file1)
outfile1 = re.sub('.fastq','', stem1) + '_withGArep.fastq.gz'
print('')
print(outfile1)
print('')
    
stem2 = re.sub('.gz', '', args.file2)
outfile2 = re.sub('.fastq','', stem2) + '_comp_GArep.fastq.gz'
print(outfile2)

# Read first file, find reads matching pattern, store reads and line 
# numbers, write first output file (GA-repeat containing).
with gzip.open(args.file1, 'rt') as handle:
    lnum = 0
    hits = []
    lnums = set()
    pattern= re.compile(args.pattern)
    for record in SeqIO.parse(handle, "fastq"):
        if pattern.search(str(record.seq)):
            hits.append(record)
            lnums.add(lnum)
        lnum += 1
        if (lnum % 1_000_000) == 0:
            print(lnum / 1_000_000, end=' ')

with gzip.open(outfile1, 'wt') as outfile:
    SeqIO.write(hits, outfile, 'fastq')

# Read second file, write reads corresponding to hits
# above to output file 2.
with gzip.open(args.file2, 'rt') as handle:
    lnum = 0
    matched_seqs = []
    for record in SeqIO.parse(handle, "fastq"):
        if lnum in lnums:
            matched_seqs.append(record)
        lnum += 1

with gzip.open(outfile2, 'wt') as outfile:
    SeqIO.write(matched_seqs, outfile, 'fastq')
"""
"""  
    