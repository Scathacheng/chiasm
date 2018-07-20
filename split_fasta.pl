use strict;
use warnings;

my $prefix = $ARGV[1];
my $max_nts = $ARGV[2];
my $head = "";
my $seq = "";
my $files = 1;
open(OUT, ">$prefix.$files.fna");
my $nts = 0;
while (<STDIN>) {
	chomp;
	if($_ =~ /^>/) {
		if($head){
			if($nts > $max_nts){
				close(OUT);
				$files++;
				open(OUT, ">$prefix.$files.fna");
				$nts = 0;
			}
			print OUT "$head\n$seq\n";
		}
		$head = $_;
		$seq = "";
	} else {
		$seq .= $_;
		$nts += length($_);
	}
}
if($head){
	print OUT "$head\n$seq\n";
}
close(OUT);
