use strict;
use warnings;

my ($a, $aa) = fasta2array($ARGV[0]);
my ($b, $bb) = fasta2array($ARGV[1]);
our @a = @$aa;
our @b = @$bb;
my $l = $ARGV[2];

open(COFOLDa,">cofoldS.$l.tmp") || die($!);
open(COFOLDb,">cofoldSE.$l.tmp") || die($!);
open(COFOLDc,">cofoldE.$l.tmp") || die($!);
open(COFOLDd,">cofoldES.$l.tmp") || die($!);
for (my $i = 0; $i < scalar(@a); $i++) {
	unless ($a[$i] && $b[$i] and length($a[$i]) > $l and length($b[$i]) > $l) {next;}
	
	my $astart = substr($a[$i],0,min(length($a[$i]),$l));
	my $bstart = substr($b[$i],0,min(length($b[$i]),$l));
	
	my $aend = reverse(substr(reverse($a[$i]),0,min(length($a[$i]),$l)));
	my $bend = reverse(substr(reverse($b[$i]),0,min(length($b[$i]),$l)));

	my $astart_cons = ( '<' x length($astart) );
	my $bstart_cons = ( '>' x length($bstart) );
	my $aend_cons = ( '<' x length($aend) );
	my $bend_cons = ( '>' x length($bend) );
	
	print COFOLDa "$astart&$bstart\n";
	print COFOLDa "$astart_cons&$bstart_cons\n";
	
	print COFOLDb "$astart&$bend\n";
	print COFOLDb "$astart_cons&$bend_cons\n";
	
	print COFOLDc "$aend&$bend\n";
	print COFOLDc "$aend_cons&$bend_cons\n";
	
	print COFOLDd "$aend&$bstart\n";
	print COFOLDd "$aend_cons&$bstart_cons\n";
}
close(COFOLDa);
close(COFOLDb);
close(COFOLDc);
close(COFOLDd);


sub min {
	if ($_[0] < $_[1]) {return $_[0];}
	return $_[1];
}

# Reads a fasta file, returns two array pointers for header an sequence
# Access e.g. via $results[x]->[y]
#my ($ph, $ps) = fasta2array($file);
#our @seqs = @$ps;
#our @header = @$ph;
sub fasta2array {
	my $file = shift;
	my $i = -1;
	my @header;
	my @sequences;
	open(FILE,"<$file") || die("fasta2array(): Could not open file: $file\n");
	while(<FILE>) {
		chomp;
		$_ =~ s/\r$//g;
		if ($_ =~ /^>/) {
			$_ =~ s/^>//;
			$i++;
			$header[$i] = $_;
			$sequences[$i] = "";
		}
		else {
			$sequences[$i] .= $_;
		}	
	}
	close(FILE);

	return(\@header,\@sequences);
}

