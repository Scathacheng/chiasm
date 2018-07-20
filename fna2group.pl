use strict;
use warnings;

my $name = $ARGV[0];
$name =~ s/mock.\d+.//;
$name =~ s/\..*//;


my $id = 0;
open(IN, "<$ARGV[0]");
while (<IN>) {
	if ($_ =~ /^>(.+?)\s/) {
		$id = $1;
		$_ =~ s/>/>$name--/;
		$id =~ s/\|/__/g;
	}
	open(OUT,">>group.$id.fasta") or die($!);
	print OUT $_;
	close(OUT);
}
close(IN);
