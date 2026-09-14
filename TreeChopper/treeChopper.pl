#!/usr/bin/env perl

use strict;
use warnings;
## disable perl's warning mechanism
no warnings 'recursion';

use FindBin;

my $usage = "usage: $0 TreeFile maxLeafDistance [link_score_cutoff]\n\n";

my $treeFile = $ARGV[0] or die $usage;
my $maxLeafDistance = $ARGV[1] or die $usage;
my $link_score = $ARGV[2];


main: {
	
	my $utilDir = "$FindBin::Bin/util/";

	## get list of all leaf nodes:
	my $leaf_file = "/tmp/tmp.$$.leafs";
	
	my $cmd = "$utilDir/tree_report_nodes.pl $treeFile > $leaf_file";
	&process_cmd($cmd);
	
	my @leaves = `cat $leaf_file`;
	unlink($leaf_file);
	chomp @leaves;
		
	## perform leaf clustering:
		
	$cmd = "$utilDir/tree_leaf_pairs_within_dist.pl $treeFile $maxLeafDistance | $utilDir/print.pl 0 1 | $utilDir/slclust";
	
	if ($link_score) {
		$cmd .= " -j $link_score";
	}
	
	my @clusters = `$cmd`;
	chomp @clusters;
	
	my %seen;
	foreach my $cluster (@clusters) {
		my @eles = split (/\s+/, $cluster);
		my $num_eles = scalar(@eles);
		print "$num_eles\t$cluster\n";
		foreach my $ele (@eles) {
			$seen{$ele} = 1;
		}
	}

	## report the singletons
	foreach my $leaf (@leaves) {
		unless ($seen{$leaf}) {
			print "1\t$leaf\n";
		}
	}


	exit(0);
}


####
sub process_cmd {
	my ($cmd) = @_;
	
	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}

	
	
