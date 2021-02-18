# Make random single-end reads from restriction digest of genome. Currently must be palindromic.
# Will need to update if using non-palindromic cutter.

from optparse import OptionParser
from Bio import SeqIO
import re
import gzip
import numpy as np

def parse_options():
	parser = OptionParser()
	parser.add_option("-g", "--genome_file", dest="genome_file",
					  help="Genome file, fasta format", metavar="GENOMEFILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="Path for output file", metavar="OUTFILE")
	parser.add_option("-n", "--reads_per_cut", dest="reads_per_cut", default=200,
					  help="Number of reads for each restriction site (cut)", metavar="N")
	parser.add_option("-c", "--cut_seq", dest="cut_seq", default='AATT',
					  help="Palindromic restriction site (default AATT)", metavar="Cutseq")
	
	(options, args) = parser.parse_args()
	return options          	

options = parse_options()
genome_file = options.genome_file
outfile = gzip.open(options.outfile + '.gz','wb')
reads_per_cut = int(options.reads_per_cut)
restriction_site = options.cut_seq

def revcomp(input):
    output = ''
    for letter in input:
        letter = letter.upper()

        if letter == 'A':
            output += 'T'
        elif letter == 'T':
            output += 'A'
        elif letter == 'G':
            output += 'C'
        else:
            output += 'G'

    return(output[::-1])

seq_id = 1

for record in SeqIO.parse(genome_file, "fasta"):
    chr_ = re.sub('chr', '', record.id)
    print(chr_)
    chr_seq = str(record.seq).upper()
    cut_sites = [m.start() for m in re.finditer(restriction_site, chr_seq)]
    for site in cut_sites:
        # Draw 5 just because.
        for n in range(0, reads_per_cut):
            shear_length = np.max([np.random.poisson(150,1)[0], 100])
            direction = np.random.choice([-1,1])
            if (direction == -1):
                start = site - shear_length
                end = site
                seq = chr_seq[start:end]
            elif (direction == 1):
                start = site
                end = start + shear_length
                seq = revcomp(chr_seq[start:end])
            if (len(seq) >= 100):
                seq = seq[:100]
                quals = '~' * 100
                line_towrite = '@' + str(seq_id) + '\n' + seq + '\n+\n' +  quals + '\n'
                outfile.write(line_towrite.encode())

            seq_id += 1
            
outfile.close()