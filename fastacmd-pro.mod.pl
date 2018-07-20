#!/usr/bin/perl
###	Fastcmd-pro
## 	Version:	1.1
##	Last updated:	2010-06-21
##	Author:		Marcus Lechner
###

# Cut multiple entries from gff or blast data
# large speedup compared to multiple fastacmd commands
# no database needed, no header restrictions (like no |, no numbers), but its picky about line lengths! (must be equal and not too long)

use Bio::DB::Fasta;
use warnings "all";
use strict;

unless (defined($ARGV[1])) {
	print STDERR "Fastacmd-pro\nExtracts sequence data form fasta based on gff or blast-m8 coordinates\nUsage: fastacmd-pro [-type=gff|blast] FASTA GFF/BLAm8 >SEQUENCES\n";
	exit;
}

my $type = "gff";

if ($ARGV[0] =~ /^-type=/) {
	$type = shift(@ARGV);
	$type =~ s/.+=//g;
	if ($type ne "gff" && $type ne "blast") {die("Error: $type is not allowed. Please use either 'gff' or 'blast'.\n");}
}


unless (-e $ARGV[0]) {
	print STDERR "Error: File '$ARGV[0]' not found!\n";
	exit;
}

unless (-e $ARGV[1]) {
	print STDERR "Error: File '$ARGV[1]' not found!\n";
	exit;
}

# Read db
my $db = Bio::DB::Fasta->new($ARGV[0],-reindex => 1);

if ($type eq "gff") {
open(GFF,"<$ARGV[1]") || die ("Error: Could not open file '$ARGV[1]'\n");
while (<GFF>) {
	if ($_ =~ /^#/) {next;}
	chomp;
	my @data = split(/\t+/,$_,8);
	unless (defined($data[6])) {
		next;
	}

	my $id = $data[0];
	my $start;
	my $stop;

	if ($data[6] eq "+") {
		$start = $data[3];
		$stop = $data[4];
	}
	else {
		$start = $data[4];
		$stop = $data[3];
	}
	
	my $should_length = abs($data[3]-$data[4])+1;

	my $seq = $db->subseq($id,$start,$stop);
	if (!defined($seq)) {
		print STDERR "Warning: $data[0] |$data[3]..$data[4]|$data[6] not found!\n";
	}
	else {
		if ($should_length != length($seq)) {print STDERR "Warning: $data[0] |$data[3]..$data[4]|$data[6] exceeds entry limits\n";next;}
		print ">$data[0]--$data[2]|$data[3]..$data[4]|$data[6]\n";
		print $seq;
		print "\n";
	}

}
close(GFF);
}
else {
open(BLA,"<$ARGV[1]") || die ("Error: Could not open file '$ARGV[1]'\n");
while (<BLA>) {
	chomp;
	my ($query_id,$subject_id,$identity,$alignment_length,$mismatches,$openings,$query_start,$query_end,$subject_start,$subject_end,$evalue,$bitscore) = split(/\s+/,$_);
	my $should_length = abs($subject_start-$subject_end)+1;
		my $header = "$query_id -> $subject_id|$subject_start..$subject_end|id:$identity,evalue:$evalue,bitscore:$bitscore";
	my $seq = $db->subseq($subject_id,$subject_start,$subject_end);
	if (!defined($seq)) {
		$subject_id =~ s/^M//;	# tmp workaround
		$header = "$query_id -> $subject_id|$subject_start..$subject_end|id:$identity,evalue:$evalue,bitscore:$bitscore";
		$seq = $db->subseq($subject_id,$subject_start,$subject_end);
		if (!defined($seq)) {
			print STDERR "Warning: $header not found!\n";
		}
	}
	if (defined($seq)) {
		if ($should_length != length($seq)) {print STDERR "Warning: $header exceeds entry limits\n";next;}
		print ">$header\n";
		print $seq;
		print "\n";
	}
}
close(BLA);
}


