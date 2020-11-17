'''
Parses a paired mapping file where the chromosome names have the species name ('mel' or 'pse') appended with an underscore.
Sorts into three files. One contains read pairs both from mel, one where both are from pseudo, and one for mixed read pairs.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Pair mapping file", metavar="FILE")

	(options, args) = parser.parse_args()
	return options
            
options = parse_options()

infile_name = options.filename
infile = open(infile_name, 'r')

file_stem = re.sub('.txt', '', infile_name)

mel_mel = open(file_stem + '_melmel.txt', 'w')
vir_vir = open(file_stem + '_virvir.txt', 'w')
mel_vir = open(file_stem + '_melvir.txt', 'w')

for line in infile:
	items = line.split()
	chr1 = items[2]
	chr2 = items[5]
	
	if (not(re.match('scaffold', chr1)) and not(re.match('scaffold', chr2)) ):
		mel_mel.write(line)
	elif (re.match('scaffold', chr1) and re.match('scaffold', chr2)):
		vir_vir.write(line)
	else:
		mel_vir.write(line)

mel_mel.close()
vir_vir.close()
mel_vir.close()
