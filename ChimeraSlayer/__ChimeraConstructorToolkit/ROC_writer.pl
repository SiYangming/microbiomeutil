#!/usr/bin/env perl

use strict;
use warnings;

## from 0-100% BS support, step 5, compute TP and FP rates from chimera and control CPC files.

my $usage = "usage: $0 chimeras.CPC  control.CPC\n\n";

my $chimera_cpc = $ARGV[0] or die $usage;
my $control_cpc = $ARGV[1] or die $usage;

print "#BS\tTP\tFP\n";
for (my $min_BS = 0; $min_BS <= 100; $min_BS += 5) {
	
	my $TP = &parse_CPC($chimera_cpc, $min_BS);
	my $FP = &parse_CPC($control_cpc, $min_BS);
	
	print join("\t", $min_BS, $TP, $FP) . "\n";
}


exit(0);


####
sub parse_CPC {
	my ($file, $min_BS) = @_;

	my $total = 0;
	my $yes = 0;

	open (my $fh, $file) or die $!;
	while (<$fh>) {
		my @x = split(/\t/);
		
		my $BS_A = $x[6];
		my $BS_B = $x[9];
		
		$total++;

		if ($BS_A >= $min_BS || $BS_B >= $min_BS) {
			$yes++;
		}
	}
	close $fh;


	return($yes/$total * 100);
}
