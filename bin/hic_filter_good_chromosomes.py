'''
Describe program here.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Pairwise mapping file", metavar="FILE")
	parser.add_option("-s", "--species", dest="species", default="mel",
					  help="Species", metavar="SPECIES")

	(options, args) = parser.parse_args()
	return options
            
options = parse_options()
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
good_chromosomes = chr_lists[species]

file = open(options.filename, 'r')

for line in file:
	line = line.rstrip()
	items = line.split()
	chr1, chr2 = items[2], items[5]
	if ((chr1 in good_chromosomes) and (chr2 in good_chromosomes)):
		print(line)
file.close()
