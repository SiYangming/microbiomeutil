#!/usr/bin/env perl

use strict;
use warnings;

open (my $fh, "rRNA16S.gold.taxonomy") or die $!;
while (<$fh>) {
	chomp;
	my ($acc, $species, $taxonomy) = split(/\t/);

	my ($genus_name, $rest) = split(/\s+/, $species, 2);
	
	my @taxons = split(/;\s+/, $taxonomy);
	
	my $pred_genus = pop @taxons;
	
	$pred_genus =~ /^(\S+)/;
	$pred_genus = $1;
	
	if ($genus_name eq $pred_genus) { 
		print "OK\t$acc\n";
	}
	else {
		print "Conflict\t$acc\t$species\t$taxonomy\n";
	}
}


exit(0);

