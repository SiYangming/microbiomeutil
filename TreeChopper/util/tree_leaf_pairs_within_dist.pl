#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'recursion';

use Bio::TreeIO;
use Data::Dumper;

my $usage = "usage: $0 treeFile max_leaf_pair_distance\n\n";

my $treeFile = $ARGV[0] or die $usage;
my $max_leaf_pair_dist = $ARGV[1] or die $usage;


main: {
	my $treeio = new Bio::TreeIO('-format' => 'newick',
								 '-file' => $treeFile);

	my $tree = $treeio->next_tree();

	my $root_node = $tree->get_root_node();
	
	my @nodes = $tree->get_nodes();
	
	## build a map of the nodes to accessions.
	my %leaf_id_to_node;
	foreach my $node (@nodes) {
		if ($node->is_Leaf()) {
			my $id = $node->id();
			$leaf_id_to_node{$id} = $node;
		}
	}
	
	
	## start from one of the nodes, and climb to the ancestor that contains it and the others as descendants.
	
	foreach my $leaf_node (values %leaf_id_to_node) {

		my %seen;
		&climb_from_node_print_leaves_within_distance($leaf_node, $leaf_node, \%seen, 0);
		
	}

	exit(0);
	
}
	

####
sub climb_from_node_print_leaves_within_distance {
	my ($climb_node, $leaf_node, $seen_href, $accumulated_distance) = @_;
	
	$seen_href->{$climb_node} = 1;
		
	$accumulated_distance += $climb_node->branch_length() || 0;
	
	if ($accumulated_distance > $max_leaf_pair_dist) {
		return;
	}
	
	## get ancestor
	my $ancestor = $climb_node->ancestor();
	if (! defined $ancestor) {
		## current node is root node.
		## cannot climb up and traverse down.
		return;
	}
	

	## walk through the descendants
	&traverse_downward($ancestor, $leaf_node, $seen_href, $accumulated_distance);
	
	## walk up one node and start down again
	&climb_from_node_print_leaves_within_distance($ancestor, $leaf_node, $seen_href, $accumulated_distance);
	return;
}


####
sub traverse_downward {
	my ($node, $leaf_node, $seen_href, $accumulated_distance) = @_;
	
	if ($node->is_Leaf()) {
		print join ("\t", $leaf_node->id(), $node->id(), $accumulated_distance) . "\n";
		return;
	}
	else {
		my @descendents = $node->each_Descendent();
		foreach my $descendent (@descendents) {
			if ($seen_href->{$descendent}) { next; }
			my $branch_length = $descendent->branch_length() || 0;
			if ($accumulated_distance + $branch_length < $max_leaf_pair_dist) {
				
				&traverse_downward($descendent, $leaf_node, $seen_href, $accumulated_distance + $branch_length);
			}
		}
		return;
	}
}
		
		


__END__

	
	my @descendents = $node->get_all_Descendents();
my %accs;
foreach my $descendent (@descendents) {
	if ($descendent->is_Leaf()) {
		my $id = $descendent->id();
		$accs{$id} = 1;
	}
}

		my $found_all_flag = 1;
		foreach my $acc (@accessions) {
			if (! exists ($accs{$acc})) {
				$found_all_flag = 0; 
				last;
			}
		}

		if ($found_all_flag) {
			print join ("\n", sort keys %accs) . "\n";
			exit(0);
		}

		$node = $node->ancestor();
	}

	## if got here, then never found a clade that contained all of them.
	
	print STDERR "Error, no clade was found to contain all accessions: @accessions\n\n"
		. "Instead, all nodes in each accessions clade from the root are provided:\n\n";

	


	#################################################################
	## Retrieving all clades containing the acessions down from the root.
	##################################################################


	my %accs = ();

	foreach my $acc (@accessions) {
		
		my $node = $acc_to_node{$acc};
		
		## climb to node connected to the root:
		while ( (my $ancestor_node = $node->ancestor()) ne $root_node) {
			$node = $ancestor_node;
		}

		my @descendants = $node->get_all_Descendents();
		foreach my $othernode ($node, @descendants) {
			if ($othernode->is_Leaf()) {
				my $id = $othernode->id();
				$accs{$id} = 1;
			}
		}

	}

	print join ("\n", sort keys %accs) . "\n";
	exit(0);
		
}


