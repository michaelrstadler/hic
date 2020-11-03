'''
Simple utility script. Takes a fasta file (presumably genome) and also a file with coordinates in UCSC format chr:pos1-pos2
and extracts the sequences, printing each one as a separate entry in a fasta file. Also can add 5' adn 3' flanking regions
v2 uses tab delmineted format
'''
from optparse import OptionParser
import sys
import re
from Bio import SeqIO

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--fasta_file", dest="fastafile",
					  help="fasta file", metavar="FILE")
	parser.add_option("-l", "--listfile", dest="listfile",
					  help="file with coordinates", metavar="LIST")
	parser.add_option("-3", "--three", dest="three",
					  help="3' extension in bp", metavar="THREE", default=0)
	parser.add_option("-5", "--five", dest="five",
					  help="5' extension in bp", metavar="FIVE", default=0)
	parser.add_option("-o", "--out_stem", dest="out_stem",
					  help="outfile stem", metavar="STEM")
	parser.add_option("-b", "--block_format", dest="block_format",
					  help="block format (no >)", metavar="BLOCK", default=False)

	(options, args) = parser.parse_args()
	return options

options = parse_options()

genome = {}
for record in SeqIO.parse(options.fastafile, "fasta"):
	genome[record.id] = record.seq

listfile = open(options.listfile, 'r')
outfile = open(options.out_stem + '_seqs_5p' + str(options.five) + '_3p' + str(options.three) + '.fa','w')
for line in listfile:
	line = line.rstrip()
	if(line[0:5] != 'track'):
		(chr, posL, posR, name, Bscore) = line.split()[0:5]
		posL = int(posL)
		posR = int(posR)
		posL = posL - int(options.five)
		posR = posR + int(options.three)
		chr_short = re.sub('chr', '', chr)

		seq = genome[chr_short][posL:posR]
		if (not options.block_format):
			outfile.write('>' + chr + '_' + str(posL) + '_' + str(posR) + '_' + Bscore + '\n')
		outfile.write(str(seq) + '\n')

listfile.close()
