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
    parser.add_option("-f", "--compressed_bin_folder", dest="compressed_bin_folder",
                      help="Compressed bin file folder", metavar="BFOLD")
    parser.add_option("-p", "--compressed_file_stem", dest="compressed_file_stem",
                      help="Compressed bin file stem", metavar="CFSTEM")
    parser.add_option("-s", "--synteny_file", dest="synteny_file",
                      help="Synteny file", metavar="SFILE")
    parser.add_option("-o", "--out_folder", dest="out_folder",
                      help="Folder for compressed bin output files", metavar="OUTFOLDER")
    parser.add_option("-r", "--ref_species", dest="ref_species",
                      help="Reference species designation within synteny file", metavar="REF")
    parser.add_option("-t", "--target_species", dest="target_species",
                      help="Target species (matches compressed bin data) designation within synteny file", 
                      metavar="TARGET")
    parser.add_option("-v", "--viewer_folder", dest="viewer_folder",
                      help="Folder for viewer files", metavar="VIEWFOLDER")
    parser.add_option("-k", "--skip_matching", dest="skip_matching",
                      help="If matched compressed bin files already made, just generate viewer files", 
                      default=False, metavar="SKIPMATCH")
    parser.add_option("-w", "--width", dest="width",
                      help="Optional: width in bins of compressed bin data", default=3500,
                      metavar="WIDTH")
    parser.add_option("-b", "--bin_size", dest="bin_size",
                      help="Optional: size of bins in bp", default=500,
                      metavar="BINSIZE")
    

    (options, args) = parser.parse_args()
    return options

match_syntenic_blocks_path = '/Users/michaelstadler/Bioinformatics/Projects/insulators/bin/match_syntenic_blocks_compressed.py'
make_viewer_files_path = '/Users/michaelstadler/Bioinformatics/Projects/insulators/bin/hic_make_viewerfiles.py'

options = parse_options()
compressed_bin_folder = options.compressed_bin_folder
synteny_file = options.synteny_file
out_folder = options.out_folder
ref_species = options.ref_species
target_species = options.target_species
width = options.width
bin_size = options.bin_size
compressed_file_stem = options.compressed_file_stem
viewer_folder = options.viewer_folder
skip_matching = options.skip_matching

files = os.listdir(compressed_bin_folder)

if (not skip_matching):
  print('woooooo')
  for f in files:
    # Avoid .DS_Store
    if (str(f)[0] != '.'):
      chromosome = re.sub(compressed_file_stem, '', f)
      chromosome = re.sub('.txt.gz', '', chromosome)
      # call block matcher
      check_call(['python', match_syntenic_blocks_path,
        '-f', os.path.join(compressed_bin_folder, f),
        '-c', chromosome,
        '-s', synteny_file,
        '-o', out_folder,
        '-r', ref_species,
        '-t', target_species,
        '-w', str(width),
        '-b', str(bin_size)])

files = os.listdir(out_folder)

for f in files:
  # Avoid .DS_Store
  if (str(f)[0] != '.'):
    chromosome = re.sub('.txt.gz', '', f)
    check_call(['python', 
      make_viewer_files_path, 
      '-f', os.path.join(out_folder, f),
      '-o', viewer_folder,
      '-c', chromosome])



