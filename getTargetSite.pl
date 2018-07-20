use strict;
use warnings;

our $ref = "amel";
our $target = $ARGV[0];
my %algn = fasta2hash($target);
if (!$algn{$ref}) {die("AMEL not part of alignment in $target\n");}
my $length = 0;
my @circ = bed2hash($ARGV[1]);
my %categories = (
	"apis" => ['ador', 'aflo', 'aces',],
	"eusocial" => ['edil', 'bter', 'bimp', 'lven', 'mqua',],
	"fly" => ['dmel',],
	"silk" => ['bmor',],
);
my %symbol = (apis => 'a', eusocial => 'e', fly => 'f', silk => 's');

my %conserved;
my %perfect;
foreach (keys %algn) {
	my @cols = split(//,$algn{$_});
	$algn{$_} = \@cols;
	if ($ref eq $_) {$length = scalar(@cols);}
	$conserved{$_} = 1;
	$perfect{$_} = 1;
}

# Find ref position
my $add_to_pos = 0;
for (my $i = 0; $i < $length; $i++) {
	if ($algn{$ref}[$i] eq "-" && $i - $add_to_pos < 100) {$add_to_pos++;}
}

# Fetch from position
for (my $i = 100+$add_to_pos; $i < 100+$add_to_pos+6; $i++) {
	foreach my $id (sort keys %algn) {
		if ($id eq $ref) {next;}
		if (uc $algn{$id}[$i] ne uc $algn{$ref}[$i]) {$perfect{$id} = 0;}
	}
}

$target =~ /--(.*)--(.*)\|(\d*)..(\d*)\|([-\+]).*/ or die('pattern not found');
my $text = "";
my $pass = 0;
foreach my $cat (sort keys %categories) {
	my $c = 0; my $p = 0;
	foreach my $g (@{$categories{$cat}}) {
		if ($conserved{$g}) {$c++;}
		if ($perfect{$g}) {$p++;}
	}
	my $max = scalar(@{$categories{$cat}});
	if($max > 1) {$max--;}
	if ($max <= $c) {
		if ($max <= $p) {$text .= uc $symbol{$cat};}
		else {$text .= $symbol{$cat};}
		$pass = 1;
	} else {$text .= "*";}
	$text .=  "\t";
}
my $n = keys %conserved;
my $start = $3+100;
my $end = $4-100;
if ($pass) {print "$1\t$start\t$end\t$2\t$n\t$5\t$text";} else {die("no conservation in $target\n");}
my @circs;
for (my $c = 0; $c < scalar(@circ); $c++) {
	my @col = @{$circ[$c]};
	if($col[0] eq $1 and $col[1] le $start+1 and $end-1 le $col[2]) {push @circs, $col[3];}
}
#if( not @circs) {die("$target does not belong to circRNA");}
print join ',', @circs;
print "\n";
exit();

open(OUT, ">$ARGV[0].target.fna");
# Fetch from position
foreach my $id (sort keys %algn) {
	print OUT ">$id\n";
	for (my $i = 100+$add_to_pos; $i < 100+$add_to_pos+6; $i++) {
		print OUT $algn{$id}[$i];
	}
	print OUT "\n";
}
close(OUT);


sub fasta2hash {
	my $file = shift;
	my $cut = shift;
	if (!defined($cut)) {$cut = 0;}
	my %hash;
	open(FILE,"<$file") || die("fasta2hash(): Could not open file: $file\n");
	my $last_hash = "";
	while(<FILE>) {
		$_ =~ s/\r/\n/g;
		chomp;
		if ($_ =~ /^>/) {
			$_ =~ s/^>//;
			if ($cut) {
				$_ =~ s/ .*$//s;
			}
			
			$_ =~ /^(....)/;
			if($1 eq $ref) {$target = $_;}
			$hash{$1} = "";
			$last_hash = $1;
			next;
		}
		$hash{$last_hash} .= $_;
	}
	close(FILE);
	return %hash;
}

sub bed2hash {
	my $file = shift;
	my @arr;
	open(FILE,"<$file") || die("gff2hash(): Could not open file: $file\n");
	my $last_hash = "";
	while(<FILE>) {
		chomp;
		if ($_ =~ /^#/ or not $_) {next;}
		push @arr, [split(/\t/)];
	}
	close(FILE);
	return @arr;
}
