# This converts my paired map file format to the format demanded for juicebox.
# it also must be sorted. I used unix sort, I think it's sort -k3,3 -k6,6 file >outfile
# I may have columns wrong, but you want to sort by chromosome 1 then chromosome 2. All incidences
# of a chromosome pair must be in one place in the file. i.e., you can't have some 2R ... 2L in one
# place and 2L .. 2R in another. This script does the work of "ordering" the chromosomes so that 
# the "lower" choromosome (arbitrary, of course) is always first. Subsequently sorting on the chromosome
# columns separately does the trick.

my %chromosome_vals = (
		"2L" => 1,
		"2R" => 2,
		"3L" => 3,
		"3R" => 4,
		"X" => 5,
		"4" => 6,
);

while(<>){
	chomp;
	my ($name, $str1, $chr1, $pos1, $str2, $chr2, $pos2) = split/\t/;


	if ($chromosome_vals{$chr1} && $chromosome_vals{$chr2}){
		$str1 = convert_strand($str1);
		$str2 = convert_strand($str2);
		($chr1, $chr2) = order_chromosomes($chr1, $chr2);
		$chr1 = 'arm_' . $chr1;
		$chr2 = 'arm_' . $chr2;
		print "$name\t$str1\t$chr1\t$pos1\t1\t$str2\t$chr2\t$pos2\t100\t1\t1\n";
	}
}

sub convert_strand{
	if ($_ eq '+'){
		return 0;
	}
	else{
		return 16;
	}
}

sub order_chromosomes{
	my ($temp1, $temp2) = ($_[0], $_[1]);
	if ($chromosome_vals{$temp1} > $chromosome_vals{$temp2}){
		return($temp2, $temp1);
	}
	else{
		return($temp1, $temp2);
	}
	
}