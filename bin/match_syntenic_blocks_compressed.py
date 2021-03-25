'''

'''
############################################################################
from optparse import OptionParser
import sys
import re
import numpy as np
import os
import sys
import gzip
from subprocess import check_call


def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--compressed_bin_file", dest="compressed_bin_file",
                      help="Compressed bin file", metavar="BFILE")
    parser.add_option("-c", "--chromosome", dest="chromosome",
                      help="Chromosome of compressed bin file", metavar="CHR")
    parser.add_option("-s", "--synteny_file", dest="synteny_file",
                      help="Synteny file", metavar="SFILE")
    parser.add_option("-o", "--out_folder", dest="out_folder",
                      help="Folder output files", metavar="OUTFOLDER")
    parser.add_option("-r", "--ref_species", dest="ref_species",
                      help="Reference species designation within synteny file", metavar="REF")
    parser.add_option("-t", "--target_species", dest="target_species",
                      help="Target species (matches compressed bin data) designation within synteny file", 
                      metavar="TARGET")
    parser.add_option("-w", "--width", dest="width",
                      help="Optional: width in bins of compressed bin data", default=3500,
                      metavar="TARGET")
    parser.add_option("-b", "--bin_size", dest="bin_size",
                      help="Optional: size of bins in bp", default=500,
                      metavar="BINSIZE")

    (options, args) = parser.parse_args()
    return options

def parse_synteny_file(file, ref_sp, target_sp):
    target = []
    ref = []
    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if (len(line) == 0):
                pass
            elif (line[0] == '#'):
                pass
            elif(line[0] == '>'):
                pass
            else:
                (species, remain) = line.split('.')
                (remain, strand) = remain.split(' ')
                (chr_, remain) = remain.split(':')
                (loc1, loc2) = remain.split('-')
                chr_ = re.sub('chr', '', chr_)
                if (species == ref_sp):
                    ref.append(chr_ + '-' + loc1 + '-' + loc2)
                elif(species == target_sp):
                    target.append([chr_, loc1, loc2, strand])
    return ref, target

def read_compressed_bin(file, width):
    init_max_size = int(1e6)
    bin_counts = np.zeros(init_max_size)
    bin_bin_counts = np.zeros((init_max_size, init_max_size))
    max_bin = 0
    if (file[-2:] == 'gz'):
        f = gzip.open(file, 'rt')
    else:
        f = open(file, 'r')

    for line in f:
        line = line.rstrip()
        if (line[0] == '#'):
            chr_, bin_, count = line[1:].split('\t')
            bin_ = int(bin_)
            count = float(count)
            bin_counts[bin_] = count
            if (bin_ > max_bin):
                max_bin = bin_
        else:
            chr_, bin1, bin2, count = line.split('\t')
            bin1 = int(bin1)
            bin2 = int(bin2)
            count = float(count)
            bin_bin_counts[bin1][bin2] = count
    f.close()
    print('Done reading' + file)
    return bin_counts, bin_bin_counts, max_bin

options = parse_options()
compressed_bin_file = options.compressed_bin_file
synteny_file = options.synteny_file
out_folder = options.out_folder
ref_species = options.ref_species
target_species = options.target_species
width = int(options.width)
chromosome = options.chromosome
bin_size = int(options.bin_size)

synteny_ref, synteny_target = parse_synteny_file(synteny_file, ref_species, target_species)
bin_counts, bin_bin_counts, max_bin = read_compressed_bin(compressed_bin_file, width)

for i in range(0, len(synteny_ref)):
    ref_name = synteny_ref[i]
    chr_, loc1, loc2, strand = synteny_target[i]
    if (chr_ == chromosome):
        bin1 = int(int(loc1) / bin_size)
        bin2 = int(int(loc2) / bin_size)
        outfile = os.path.join(out_folder, ref_name + '.txt')
        with open(outfile, 'w') as out:
            sub_bin_counts = bin_counts[bin1:bin2]
            sub_bin_bin_counts = bin_bin_counts[bin1:bin2, bin1:bin2]
            # Flip if negative strand
            if (strand == '-'):
                sub_bin_counts = sub_bin_counts[::-1]
                sub_bin_bin_counts = sub_bin_bin_counts[::-1, ::-1]
            for j in range(0, len(sub_bin_counts)):
                count = sub_bin_counts[j]
                if (count > 0):
                    out.write('#' + chromosome + '\t' + str(j) + '\t' + str(count) + '\n')

            for j in range(0, len(sub_bin_counts)):
                left = np.max((0, j - width))
                right = np.min((j + width, len(sub_bin_counts)))
                for k in range(left, right):
                    count = sub_bin_bin_counts[j,k]
                    if (count > 0):
                        out.write(chromosome + '\t' + str(j) + '\t' + str(k) + '\t'+ str(count) + '\n')
        check_call(['gzip', outfile])
            