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
pse_pse = open(file_stem + '_psepse.txt', 'w')
mel_pse = open(file_stem + '_melpse.txt', 'w')

for line in infile:
	items = line.split()
	print(items)
	chr1 = items[2]
	chr2 = items[5]
	species1 = chr1[-3:]
	species2 = chr2[-3:]
	if (species1 == 'mel' and species2 == 'mel'):
		mel_mel.write(line)
	elif (species1 == 'pse' and species2 == 'pse'):
		pse_pse.write(line)
	elif ((species1 == 'pse' and species2 == 'mel') or (species1 == 'mel' and species2 == 'pse')):
		mel_pse.write(line)

mel_mel.close()
pse_pse.close()
mel_pse.close()
