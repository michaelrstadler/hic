# -*- coding: utf-8 -*-
"""
This version ultimately failed, but it's interesting enough to save. Basically it reads parallel chunks
of the two files. As long as pairs are within the chunk size of each other, they are guaranteed
to be successfully paired. Unfortunately in practice I found that pairs could sometimes be separated 
by quite large distances and so I would not recover all pairs successfully.

@author: MStadler
"""
import sys
from optparse import OptionParser
import re
import gzip
from time import sleep
import gc

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--file1", dest="filename1",
                      help="First mapping file", metavar="FILE1")
    parser.add_option("-s", "--file2", dest="filename2",
                      help="Second mapping file", metavar="FILE2")
    parser.add_option("-c", "--chunksize", dest="chunk_size", default = 50_000_000,
                      help="Size in lines of file chunks to read", metavar="CHUNKSIZE")

    (options, args) = parser.parse_args()
    return options

def Get5prime(splitline):
    """Find 5' position of read."""
    strand = splitline[1]
    chr = splitline[2]
    pos = splitline[3]
    seq = splitline[4]
    if (strand == '+'):
        pass
    else:
        size = len(seq)
        pos = str(int(pos) + size - 1)
    return(pos)

def read_file_chunk(handle, n):
    data = {}
    count = 0
    for line in handle:
        count += 1
        pair_name = line.split()[0]
        splitline= line.split('\t')
        pos5p = Get5prime(splitline)
        linedata = "\t".join(splitline[1:3] + [pos5p])
        data[pair_name] = linedata
        if (count == n):
            break
    return data, count

options = parse_options()
# Initialize data structures and set outfile path.
file1_data = {}
chunk_size = int(options.chunk_size)
goodPair_count = 0
total = 0
distance_minimum = 0
outfilename = ''
if (re.search('R1', options.filename1)):
    outfilename = options.filename1
elif (re.search('R1', options.filename2)):
    outfilename = options.filename2
outfilename = re.sub('_R1_','_', outfilename)
outfilename = re.sub('.gz', '', outfilename)
outfilename = re.sub('.bowtie', '', outfilename)
outfilename = outfilename +  '_pairMerged_new.txt'
outfile = open(outfilename, "w")

def open_filehandle(filename):
    if (filename[-2:] == 'gz'):
            return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

file1 = open_filehandle(options.filename1)
file2 = open_filehandle(options.filename2)

def match(data1, data2, matched_reads, outfile):
    for readname in data1:
        if readname in data2:
            splitline1 = data1[readname].split('\t')
            splitline2 = data2[readname].split('\t')
            outfile.write(readname + "\t" + splitline1[0] + '\t' + splitline1[1] + '\t' + splitline1[2] + '\t' + splitline2[0] + '\t' + splitline2[1] + '\t' + splitline2[2] + '\n')
            matched_reads.append(readname)
    
data_file1 = {}
data_file2 = {}
good_pair_count = 0
total_count = 0
data_oldchunk_file1 = {}
data_oldchunk_file2 = {}
chunksize_file1 = 1
chunksize_file2 = 1
while (chunksize_file1 > 0) or (chunksize_file2 > 0):
    # Read file1 chunk, add to data.
    data_newchunk_file1, chunksize_file1 = read_file_chunk(file1, chunk_size) 
    data_newchunk_file2, chunksize_file2 = read_file_chunk(file2, chunk_size)

    print(chunksize_file2)
    #quit()
    total_count += chunksize_file1
    matched_reads = []
    #match(data_oldchunk_file1, data_oldchunk_file2, matched_reads, outfile)
    match(data_oldchunk_file1, data_newchunk_file2, matched_reads, outfile)
    match(data_newchunk_file1, data_oldchunk_file2, matched_reads, outfile)
    match(data_newchunk_file1, data_newchunk_file2, matched_reads, outfile)
    #match(data_oldchunk_file2, data_newchunk_file1, matched_reads, outfile)
    #match(data_newchunk_file2, data_oldchunk_file1, matched_reads, outfile)
    #match(data_newchunk_file2, data_newchunk_file1, matched_reads, outfile)
    
    good_pair_count += len(matched_reads)     
    for data in [data_newchunk_file1, data_newchunk_file2]: 
        for readname in matched_reads:
            if readname in data:
                del data[readname]
    #print(len(data_file1))
    del data_oldchunk_file1
    del data_oldchunk_file2
    data_oldchunk_file1 = data_newchunk_file1
    data_oldchunk_file2 = data_newchunk_file2
    gc.collect()
    
pct = float(good_pair_count) / float(total_count) * 100
print(str(good_pair_count) + ' of ' + str(total_count) + ' reads have aligned pairs, or ' + str(pct) + ' percent' )

# Read file2 chunk, add to data.
# Run over all file1 entries, look for matches in file2
# For matches, print and delete entry in dict.
# Die when the size of both chunks is zero.

"""
f2 = options.filename2
if (f2[-2:] == 'gz'):
        file2 = gzip.open(f2, 'rt')
else:
    file2 = open(f2, 'r')

for line in file2:
    total += 1
    pair_name = line.split()[0]
    splitline2= line.split('\t')
    if pair_name in file1_data:
        splitline1 = file1_data[pair_name].split('\t')
        pos_5p_2 = Get5prime(splitline2)
        outfile.write(pair_name + "\t" + splitline1[0] + '\t' + splitline1[1] + '\t' + splitline1[2] + '\t' + splitline2[1] + '\t' + splitline2[2] + '\t' + pos_5p_2 + '\n')
        goodPair_count += 1
"""

file2.close()
file2.close()


#pct = float(goodPair_count) / float(total) * 100
#print >> sys.stderr, str(goodPair_count) + ' of ' + str(total) + ' reads have aligned pairs, or ' + str(pct) + ' percent'
