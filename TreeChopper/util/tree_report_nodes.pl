#!/usr/bin/env perl

use strict;
use warnings;

use Bio::TreeIO;
use Data::Dumper;

my $usage = "usage: $0 treeFile\n\n";

my $treeFile = $ARGV[0] or die $usage;

main: {
	my $treeio = new Bio::TreeIO('-format' => 'newick',
								 '-file' => $treeFile);

	my $tree = $treeio->next_tree();

	my @nodes = $tree->get_nodes();
	
	## build a map of the nodes to accessions.
	my %acc_to_node;
	foreach my $node (@nodes) {
		if ($node->is_Leaf()) {
			my $id = $node->id();
			print "$id\n";
		}
	}
	
	exit(0);
			
}


