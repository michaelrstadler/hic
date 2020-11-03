'''
number of members of second file that are within d of features in first file
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--first_file", dest="first_filename",
					  help="first file", metavar="FIRSTFILE")
	parser.add_option("-s", "--second_file", dest="second_filename",
					  help="second file", metavar="SECONDFILE")
	parser.add_option("-d", "--distance", dest="distance",
					  help="max distance to count as overlap", metavar="DIST", default = 1000)


	(options, args) = parser.parse_args()
	return options

def overlaps(chr, point, features, sign, distance):
	for i in range(0, distance):
		new_point = point + (sign * i)
		if (chr in features):
			if(new_point in features[chr]):
				return True
	return False

def fix_chr(original):
	if (original[:3] != 'chr'):
		return('chr' + original)
	return(original)
            
options = parse_options()

distance_cutoff = int(options.distance)
features_1 = {"chr2L":{}, "chr2R":{}, "chr3L":{}, "chr3R":{}, "chrX":{}, "chr4":{}}
overlap_count = 0
total_count = 0
distance = int(options.distance)

firstfile = open(options.first_filename, 'r')
for line in firstfile:
	line = line.rstrip()
	items = line.split()
	#print(items)
	midpoint = int((float(items[1]) + float(items[2])) / 2)
	chr = fix_chr(items[0])
	if (chr in features_1):
		features_1[chr][midpoint] = 1

firstfile.close()

secondfile = open (options.second_filename, 'r')
for line in secondfile:
	total_count += 1
	line = line.rstrip()
	items = line.split()
	chr = fix_chr(items[0])
	midpoint = int((float(items[1]) + float(items[2])) / 2)
	if(overlaps(chr, midpoint, features_1, -1, distance) or overlaps(chr, midpoint, features_1, 1, distance)):
		overlap_count += 1

secondfile.close()

print(str(overlap_count) + '/' + str(total_count))
