"""

"""

from optparse import OptionParser
import sys
import re
import gzip
import io
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="paired alignment files, comma separated", metavar="FILE")
	
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")

	parser.add_option("-s", "--species",
						dest="species", default='mel',
					  help="species: mel, pseudo, or vir")

	parser.add_option("-m", "--max_chr_size",
						dest="max_chr_size", default=1e8,
					  help="Max size to allocate for each chromosome. For many chromosomes (scaffolds), can get memory intensive if set high.")

	parser.add_option("-o", "--output_stem",
						dest="file_stem", default='none',
					  help="output file stem. Adds _binCounts_Xkb.txt")

	(options, args) = parser.parse_args()
	return options
                
options = parse_options()
bin_size = int(options.bin_size)
bin_counts = {}
max_bins = {}
filenames = options.filenames
files = filenames.split(',')
line_count = 0
max_chr_size = int(options.max_chr_size)
min_chr_size = 1e6
max_chr_bins = int(max_chr_size // bin_size)
min_chr_bins = int(min_chr_size // bin_size)
species = options.species
chr_lists = {
	'2L': ['2L'],
	'mel': ['2L', '2R', '3L', '3R', 'X', '4', 'Y'],
	'vir': ['scaffold_13049','scaffold_12963','scaffold_13047','scaffold_12875','scaffold_12855',
	'scaffold_12970','scaffold_12723','scaffold_12928','scaffold_12822','scaffold_13042','scaffold_12958',
	'scaffold_12726','scaffold_13246','scaffold_12823','scaffold_13324','scaffold_12932','scaffold_13052',
	'scaffold_12758','scaffold_13050','scaffold_10324','scaffold_12728','scaffold_12954','scaffold_12472',
	'scaffold_12736','scaffold_12967','scaffold_13045','scaffold_13036','scaffold_12936','scaffold_12799',
	'scaffold_10322','scaffold_12100','scaffold_12734','scaffold_12937','scaffold_12929','scaffold_12941',
	'scaffold_12528','scaffold_12735','scaffold_12815','scaffold_12942','scaffold_10419','scaffold_9476',
	'scaffold_12847','scaffold_12930','scaffold_12790','scaffold_12730','scaffold_12891','scaffold_12952',
	'scaffold_12961','scaffold_12325','scaffold_12960','scaffold_12945','scaffold_12514']
}
chromosomes = chr_lists[species]

def add_read(chr1, chr2, bin1, bin2, bin_counts, max_chr_bins, chromosomes):
	if ((chr1 in chromosomes) and (chr2 in chromosomes)):
		if (chr1 not in bin_counts):
			bin_counts[chr1] = {}
		if (chr2 not in bin_counts[chr1]):
			bin_counts[chr1][chr2] = np.zeros((max_chr_bins, max_chr_bins))
		bin_counts[chr1][chr2][bin1][bin2] += 1


def update_max(chr_, bin_, max_bins):
	if (chr_ not in max_bins):
		max_bins[chr_] = bin_
	elif(bin_ > max_bins[chr_]):
		max_bins[chr_] = bin_

for f in files:
	print('Reading ' + f + '...')
	if (f[-2:] == 'gz'):
		file1 = gzip.open(f, 'rt')
	else:
		file1 = open(f, 'r')

	for line in file1:
		line_count = line_count + 1
		if (line_count % 10000000 == 0):
			print('.')
		line = line.rstrip()
		items = line.split()
		(chr1, Lmost1, chr2, Lmost2) = items[2], int(items[3]), items[5], int(items[6])
		bin1 = int(Lmost1 / bin_size)
		bin2 = int(Lmost2 / bin_size) 
		add_read(chr1, chr2, bin1, bin2, bin_counts, max_chr_bins, chromosomes)

		add_read(chr2, chr1, bin2, bin1, bin_counts, max_chr_bins, chromosomes)
		update_max(chr1, bin1, max_bins)
		update_max(chr2, bin2, max_bins)
	file1.close()
file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

#outfile = open(file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt','w')
outfilename = file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt'


#print column labels
"""
for chr3 in chromosomes:
	for bin3 in range(0, max_bins[chr3] + 1):
		outfile.write('\tchr' + chr3 + '_' + str(bin3))
outfile.write('\n')
"""

# Save data for all chromosomes.
for chr1 in chromosomes:

	data_this_chr = bin_counts[chr1][chromosomes[0]][0:(max_bins[chr1]+1), 0:(max_bins[chromosomes[0]]+1)]
	
	for chr2 in chromosomes[1:]:
		if (chr2 in bin_counts[chr1]):
			data_this_chr = np.hstack((data_this_chr, bin_counts[chr1][chr2][0:(max_bins[chr1]+1),0:(max_bins[chr2]+1)]))
		else:
			data_this_chr = np.hstack((data_this_chr, np.zeros((max_bins[chr1]+1,max_bins[chr2]+1))))
	if (chr1 == chromosomes[0]):
		data_all = data_this_chr
	else:
		data_all = np.vstack((data_all, data_this_chr))

np.savetxt(outfilename, data_all, newline='\n', fmt='%.6e')
