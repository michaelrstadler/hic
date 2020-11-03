# Quick perly work to take a bed file and convert it to pseudo-gff file. I had some scripts
# that needed gff files and could only download BED. This was made for converting gene structure BED
# to a gff format for mapping distances from boundaries to TSS (promoters).

use warnings;
use strict;

while(<>){
	chomp;
	my @line = split/\t/;
	if ($line[5] eq '+'){
		$line[2] = $line[1] + 1;
	}
	if ($line[5] eq '-'){
		$line[1] = $line[2];
		$line[2]++
	}
	#print join("\t",@line);
		#print "\n";
	print "$line[0]\t.\t.\t$line[1]\t$line[2]\t.\t.\t.\t.\n";
}