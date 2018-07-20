use strict;
use warnings;

my %hash;
while (<STDIN>) {
	my ($q) = split;
	if ($hash{$q}) {next;}
	$hash{$q}++;
	print $_;
}
