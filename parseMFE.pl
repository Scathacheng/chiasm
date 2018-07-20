use strict;
use warnings;

my $gc = 0;
my $gu = 0;
my $n = 0;
while (<STDIN>) {
	if ($_ =~ /\(\s*(-\d+\.\d+)\)/) {
        print "$1\t$gc\t$gu\n";
        $gc = 0;
        $gu = 0;
        $n = 0;
    } else {
        chomp;
        $n =()= ($_ =~ m/(N|n)/g);
        $gc =()= ($_ =~ m/(G|C|g|c)/g);
        $gu =()= ($_ =~ m/(G|U|T|g|u|t)/g);
        if($n == length($_)-1) {$gc=0; $gu=0; next;}
        $gc /= length($_)-$n-1;
        $gu /= length($_)-$n-1;
    } 
}
