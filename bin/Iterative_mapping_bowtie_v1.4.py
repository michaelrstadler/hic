"""
Wrapper for mapping Hi-C reads. The challenge with standard mapping is that, particularly for 
high-resolution cutters, long reads are likely to run into Hi-C ligation junctions and thus
fail to map. Iterative mapping addresses this by working from the 5' end of read1 and read2,
starting with the first 20 bp, and adding more and more until a unique mapping can be established.

It begins by mapping the first n bp of the reads (default=20).
Uniquely mapping reads (--best --strata -m1 tags) are kept and written to the output file,
unmapped reads are writting to an unmapped reads file, and reads that map repetitively are
kept. The repetitively mapping reads and then mapped again, this time using the first n + x 
(default x=7, so 27 bp for round 2). This process is repeated until the entire read is mapped.

v1.2: works on gzipped files
v1.3: does iterative mapping up to 50 bp, then just maps entire read
"""

from optparse import OptionParser
import sys
import subprocess
from subprocess import call
import re
version = "1_4"

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
max_length_for_iterating = 50

file_stem = options.outfile_stem + '_'
num_cores = '-p' + options.cores

# Do initial round of mapping as standalone; capture unmapped reads here
file_to_map = fastq_file
unmapped_file = file_stem + 'mappingV' + version +  "_unmapped.fastq"
multiplyMapped_file = file_stem + "multiplyMapped_" + str(init_truncation_length) + ".fastq"
unique_file = file_stem + "unique_" + str(init_truncation_length) + ".bowtie"


if (file_to_map[-2:] == 'gz'):
	p1 = subprocess.Popen(["gzip", "-cd", file_to_map], stdout=subprocess.PIPE)
	call([ "./bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim),  "--un", unmapped_file, "--max", multiplyMapped_file, genome_path, "-", unique_file], stdin=p1.stdout)
else:
	call(["bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim),  "--un", unmapped_file, "--max", multiplyMapped_file, genome_path, file_to_map, unique_file])

cat_command = ["cat", unique_file]

# Build truncation lengths to use
truncation_lengths = []
for k in range(init_truncation_length + iterator, max_length_for_iterating, iterator):
	truncation_lengths.append(k)
truncation_lengths.append(read_length)

# Call bowtie for all trunctation lengths
for truncation_length in truncation_lengths:
	threeprime_trim = read_length - truncation_length
	file_to_map = multiplyMapped_file
	multiplyMapped_file = file_stem + 'mappingV' + version + "_multiplyMapped_" + str(truncation_length) + ".fastq"
	unique_file = file_stem + "unique_" + str(truncation_length) + ".bowtie"
	call(["bowtie", num_cores, "-m1", "--best", "--strata", "--trim3", str(threeprime_trim), "--max", multiplyMapped_file, genome_path, file_to_map, unique_file])
	call(["rm", file_to_map])
	cat_command.append(unique_file)
	
# Clean up files
unique_combined = open(file_stem + 'mappingV' + version + '_unique.bowtie', "w")
call(cat_command, stdout = unique_combined)
cat_command[0] = "rm"
call(cat_command)