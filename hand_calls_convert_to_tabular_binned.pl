#converts hand-called boundaries to tabular form readable by other programs. Converts to nearest 500 bp bin.

use warnings;
use strict;

my %list;

while(<>){
	chomp;
	my @thing = split/:/;
	my $pos1 = int($thing[1] / 500) * 500;
	my $pos2 = $pos1 + 499;
	my $chr = 'chr' . $thing[0];
	#print "$chr\t$pos1\t$pos2\tboundary\t1\n";
	my $name = $chr . $pos1 . $pos2;
	$list{$name} = "$chr\t$pos1\t$pos2\tboundary\t1\n";
}

foreach my $item (keys %list){
	print "$list{$item}";
}