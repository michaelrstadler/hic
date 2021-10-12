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
chr_counts = {}
filenames = options.filenames
files = filenames.split(',')
line_count = 0

species = options.species
chr_lists = {
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

def add_read(chr1, chr2, chr_counts, chromosomes):
	if ((chr1 in chromosomes) and (chr2 in chromosomes)):
		if (chr1 not in chr_counts):
			chr_counts[chr1] = {}
		if (chr2 not in chr_counts[chr1]):
			chr_counts[chr1][chr2] = 0
		chr_counts[chr1][chr2] += 1


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

		add_read(chr1, chr2, chr_counts, chromosomes)
		if chr1 != chr2:
			add_read(chr2, chr1, chr_counts, chromosomes)

	file1.close()
file_stem = ''

if (options.file_stem == 'none'):
	file_stem = re.sub('.txt', '', files[0])
else:
	file_stem = options.file_stem

#outfile = open(file_stem + '_binCounts_' + str(bin_size // 1000) + 'kB.txt','w')
outfilename = file_stem + '_chrCounts.txt'
outfile = open(outfilename, 'w')

#print column labels
outfile.write(chromosomes[0])

for chr3 in chromosomes[1:]:
	outfile.write('\t' + chr3)
outfile.write('\n')


#chromosomes = filter_chr_size(list(bin_counts.keys()), max_bins, min_chr_bins)

for chr1 in chromosomes:
	outfile.write(chr1)
	for chr2 in chromosomes:
		try:
			outfile.write('\t' + str(chr_counts[chr1][chr2]))
		except:
			outfile.write('\t0')
	outfile.write('\n')

