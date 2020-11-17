'''
Takes a list of machine-defined and hand-picked boundaries. For each machine boundary, finds distance to nearest
hand-called boundary and prints it. Hard-coded for 500 bp bins right now. Hand file needs to be binned, not raw.
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-m", "--machine_file", dest="machine_filename",
					  help="boundaries picked by algorithm", metavar="MFILE")
	parser.add_option("-f", "--hand_file", dest="hand_filename",
					  help="boundaries picked by hand", metavar="HFILE")
	parser.add_option("-d", "--distance", dest="distance",
					  help="max distance to merge boundaries", metavar="DIST", default = 1000)


	(options, args) = parser.parse_args()
	return options

def dist_to_nearest(chr, point, hand_boundaries, sign):
	for i in range(0, 1000000, 500):
		new_point = point + (sign * i)
		if(new_point in hand_boundaries[chr]):
			return (new_point - point)
	return 10000000

            
options = parse_options()

distance_cutoff = int(options.distance)
hand_boundaries = {"chr2L":{}, "chr2R":{}, "chr3L":{}, "chr3R":{}, "chrX":{}, "chr4":{}}

handfile = open(options.hand_filename, 'r')
for line in handfile:
	line = line.rstrip()
	items = line.split()
	midpoint = (float(items[1]) + float(items[2])) / 2
	chr = items[0]
	hand_boundaries[chr][midpoint] = 1

handfile.close()

machinefile = open (options.machine_filename, 'r')
for line in machinefile:
	line = line.rstrip()
	items = line.split()
	chr = items[0]
	midpoint = (float(items[1]) + float(items[2])) / 2
	dist_L = dist_to_nearest(chr, midpoint, hand_boundaries, -1)
	dist_R = dist_to_nearest(chr, midpoint, hand_boundaries, 1)
	dist = 0
	if (abs(dist_L) < abs(dist_R)):
		dist = dist_L
	else: dist = dist_R
	
	if (abs(dist) <= distance_cutoff):
		#print(items[0] + '\t' + str(items[1] + dist) + '\t' + str(items[2] + dist) + '\t' + items[3] + '\t' + str(items[4]) + '\t' + str(abs(dist)))
		new_Lmost = int(int(items[1]) + dist)
		new_Rmost = int(int(items[2]) + dist)
		print(items[0] + '\t' + str(new_Lmost) + '\t' + str(new_Rmost) + '\t' + items[3] + '\t' + str(items[4]) + '\t' + str(abs(dist)))

machinefile.close()
