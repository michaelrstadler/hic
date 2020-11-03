'''
Describe program here.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-l", "--locations", dest="loc_file",
					  help="File of locations to map", metavar="FILE2")

	(options, args) = parser.parse_args()
	return options
            
options = parse_options()
