"""
Takes a fastq file of untrimmed reads and iteratively maps, starting with short reads,
keeping reads that uniquely map, discarding those that fail to map, and then repeating
mapping of multiply-mapping reads with less truncation. This process is repeated at
intervals until the whole read is mapped.

"""

from optparse import OptionParser
import sys
import subprocess
from subprocess import call
import re
version = "1_0"
#from subprocess import check_output

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="fastq file of untrimmed reads", metavar="FILE")
	parser.add_option("-i", "--initial_length",
					   dest="initial_length", default=20,
					  help="initial truncated read length")
	parser.add_option("-l", "--read_length",
					  dest="read_length",
					  help="length of the reads supplied")
	parser.add_option("-g", "--genome",
					  dest="genome",
					  help="path of genome to map against")
	
	parser.add_option("-o", "--outfile_stem",
					  dest="outfile_stem",
					  help="stem for output files")
					  
	parser.add_option("-p", "--cores",
					  dest="cores",
					  help="number of cores to use")

	(options, args) = parser.parse_args()
	return options
	
	
	
options = parse_options()

init_truncation_length = int(options.initial_length)
read_length = int(options.read_length)
threeprime_trim = read_length - init_truncation_length
genome_path = options.genome
fastq_file = options.filename
iterator = 7

#m = re.search('.fastq', fastq_file)
file_stem = options.outfile_stem + '_'
num_cores = '-p' + options.cores

# Do initial round of mapping as standalone because we capture unmapped reads here

file_to_map = fastq_file
unmapped_file = file_stem + 'mappingV' + version +  "_unmapped.fastq"
multiplyMapped_file = file_stem + "multiplyMapped_" + str(init_truncation_length) + ".fastq"
unique_file = file_stem + "unique_" + str(init_truncation_length) + ".bowtie"


if (file_to_map[-2:] == 'gz'):
	call(["gzip", "-cd", file_to_map, "|", "bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim),  "--un", unmapped_file, "--max", multiplyMapped_file, genome_path, "-", unique_file])

else:
	call(["bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim),  "--un", unmapped_file, "--max", multiplyMapped_file, genome_path, file_to_map, unique_file])

# Successively truncate and map

truncation_length = init_truncation_length + iterator
threeprime_trim = read_length - truncation_length
cat_command = ["cat", unique_file]

while truncation_length < read_length:
	file_to_map = multiplyMapped_file
	multiplyMapped_file = file_stem + 'mappingV' + version + "_multiplyMapped_" + str(truncation_length) + ".fastq"
	unique_file = file_stem + "unique_" + str(truncation_length) + ".bowtie"
	#call(["./bowtie", "-p2", "-m1", "--best", "--strata", " --trim3 ", str(threeprime_trim), " --un ", unmapped_file,  " --max ", multiplyMapped_file, genome_path, " 2> ", stats_file, " >", unique_map_file]) 
	call(["bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim), "--max", multiplyMapped_file, genome_path, file_to_map, unique_file])
	truncation_length = truncation_length + iterator
	threeprime_trim = read_length - truncation_length
	call(["rm", file_to_map])
	cat_command.append(unique_file)
	
# Clean up files and get stats
unique_combined = open(file_stem + 'mappingV' + version + '_unique.bowtie', "w")
call(cat_command, stdout = unique_combined)
cat_command[0] = "rm"
call(cat_command)