use strict;
use warnings;

my ($a, $aa) = fasta2array($ARGV[0]);
my ($b, $bb) = fasta2array($ARGV[1]);
my @a = @$a;
my @b = @$b;
my @aa = @$aa;
my @bb = @$bb;


open(STARTS, ">blast_S.tmp");
open(ENDS, ">blast_E.tmp");
open(SE, ">blast_SE.tmp");
open(ES, ">blast_ES.tmp");
my $r = '';
for (my $i = 0; $i < scalar(@a); $i++) {
	unless ($aa[$i] && $bb[$i]) {next;}
	
	my $astart = substr($aa[$i],0,min(length($aa[$i]),100));
	my $bstart = substr($bb[$i],0,min(length($bb[$i]),100));
	
	my $aend = reverse(substr(reverse($aa[$i]),0,min(length($aa[$i]),100)));
	my $bend = reverse(substr(reverse($bb[$i]),0,min(length($bb[$i]),100)));
	
	open(STARTa, ">as.tmp");
	print STARTa ">$a[$i]\n$astart\n";
	close(STARTa);
	open(STARTb, ">bs.tmp");
	print STARTb ">$b[$i]\n$bstart\n";
	close(STARTb);
	open(ENDa, ">ae.tmp");
	print ENDa ">$a[$i]\n$aend\n";
	close(ENDa);
	open(ENDb, ">be.tmp");
	print ENDb ">$b[$i]\n$bend\n";
	close(ENDb);
	
	$r = `blastn -query as.tmp -subject bs.tmp -outfmt 6 -word_size 6 | awk '\$10 < \$9{print}' | cut -f12 | sort -r | head -n1`;
	chomp $r;
	if(not $r) {$r=0.0;}
	print STARTS "$r\n";
	
	$r = `blastn -query ae.tmp -subject be.tmp -outfmt 6 -word_size 6 | awk '\$10 < \$9{print}' | cut -f12 | sort -r | head -n1`;
	chomp $r;
	if(not $r) {$r=0.0;}
	print ENDS "$r\n";

	$r = `blastn -query as.tmp -subject be.tmp -outfmt 6 -word_size 6 | awk '\$10 < \$9{print}' | cut -f12 | sort -r | head -n1`;
	chomp $r;
	if(not $r) {$r=0.0;}
	print SE "$r\n";
	
	$r = `blastn -query ae.tmp -subject bs.tmp -outfmt 6 -word_size 6 | awk '\$10 < \$9{print}' | cut -f12 | sort -r | head -n1`;
	chomp $r;
	if(not $r) {$r=0.0;}
	print ES "$r\n";
}
close(STARTS);
close(ENDS);
close(SE);
close(ES);


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

