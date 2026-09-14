#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;

use lib ("$FindBin::Bin/PerlLib");
use Fasta_reader;
use CdbTools;
use File::Basename;

my $db_NAST = "$FindBin::Bin/../RESOURCES/rRNA16S.gold.NAST_ALIGNED.fasta";
my $db_FASTA = "$FindBin::Bin/../RESOURCES/rRNA16S.gold.fasta";


my $usage = <<_EOUSAGE_;

##########################################################################################
#
#  Required:
#
#    --query_NAST      multi-fasta file containing query sequences in alignment format
#      or
#    --query_FASTA     will first run NAST-iEr to generate NAST alignments.
#
#  Common opts:
#
#    --db_NAST        db in NAST format (default: $db_NAST)
#    --db_FASTA       db in fasta format (megablast formatted) (default: $db_FASTA)
#
#
#    -n       number of top matching database sequences to compare to (default 15)
#    -R       min divergence ratio default: 1.007
#    -P       min percent identity among matching sequences (default: 90)
#
#  ## parameters to tune ChimeraParentSelector:
#   
#  Scoring parameters:
#   -M match score   (default: +5)
#   -N mismatch penalty  (default: -4)
#   -Q min query coverage by matching database sequence (default: 70)
#   -T maximum traverses of the multiple alignment  (default: 1)

#
#  ## parameters to tune ChimeraPhyloChecker:
#
#
#    --windowSize                default 50
#    --windowStep                default 5
#    --minBS      minimum bootstrap support for calling chimera (default: 90)
#    --num_BS_replicates         default: 100
#    --low_range_finer_BS (default: 10)    If computed BS is between minBS and (minBS - low_range_finer_BS), then num_finer_BS_replicates computed.
#    --num_finer_BS_replicates (default: 1000)
#    -S           percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 10)
#    --num_parents_test       number of potential parents to test for chimeras (default: 3)
#    --MAX_CHIMERA_PARENT_PER_ID    Chimera/Parent alignments with perID above this are considered non-chimeras (default 100; turned off) 
#    
#
#  ## misc opts
#
#   --printFinalAlignments          shows alignment between query sequence and pair of candidate chimera parents
#   --printCSalignments             print ChimeraSlayer alignments in ChimeraSlayer output
#   --exec_dir                      chdir to here before running
#
#  ## testing options
#   --ignore_parent_same_acc_substring          ignores parent sequence where parent acc is a substring of the query                              
#   --ignore_identical_parent_seq               ignore parent sequences that are entirely identical to the query                                  
# 
#########################################################################################

_EOUSAGE_
	
	;

## required options
my ($query_NAST);
my $query_FASTA;

# common opts
my $numSeqsCompare = 15;
my $minPerID = 90;
my $minDivR = 1.007;

# chimera parent selector opts:
my ($matchScore, $mismatchPenalty, $maxTraverse, $minQueryCoverage);

# chimera phylochecker opts:
my $minBS = 90;
my ($windowSize, $windowStep, $minPercentSampleSNPs, $num_parents_test, $num_BS_replicates,
	$low_range_finer_BS, $num_finer_BS_replicates,
);

# misc
my $printFinalAlignments = 0;
my $help;
my $exec_dir;
my $printCSalignments;
my $MAX_CHIMERA_PARENT_PER_ID;

                                                                                                                                                
my $IGNORE_PARENT_SAME_ACC_SUBSTRING = 0;                                                                                                         
my $IGNORE_IDENTICAL_PARENT_SEQ = 0;

&GetOptions ("query_NAST=s" => \$query_NAST,
			 "query_FASTA=s" => \$query_FASTA,

			 "db_NAST=s" => \$db_NAST,
			 "db_FASTA=s" => \$db_FASTA,
			 
			 # common opts
			 "R=f" => \$minDivR,  # same as -R on chimeraMaligner
			 "n=i" => \$numSeqsCompare,
			 "P=f" => \$minPerID,
			 
			 ## ChimeraParentSelector parameters:
			 "M=i" => \$matchScore,
			 "N=i" => \$mismatchPenalty,
			 "Q=i" => \$minQueryCoverage,
			 "T=i" => \$maxTraverse,

			 ## ChimeraPhyloChecker parameters
			 			 
			 "windowSize=i" => \$windowSize,
			 "windowStep=i" => \$windowStep,
			 "minBS=i" => \$minBS,
			 'num_BS_replicates=i' => \$num_BS_replicates,
			 'low_range_finer_BS=i' => \$low_range_finer_BS,
			 'num_finer_BS_replicates=i' => \$num_finer_BS_replicates,
			 'S=i' => \$minPercentSampleSNPs,
			 "num_parents_test=i" => \$num_parents_test,
			 "MAX_CHIMERA_PARENT_PER_ID=f" => \$MAX_CHIMERA_PARENT_PER_ID,
			 
			 # misc

			 "printFinalAlignments" => \$printFinalAlignments,
			 "exec_dir=s" => \$exec_dir,
			 "printCSalignments" => \$printCSalignments,

			 # testing opts
			 'ignore_parent_same_acc_substring' => \$IGNORE_PARENT_SAME_ACC_SUBSTRING,                                                            
             'ignore_identical_parent_seq' => \$IGNORE_IDENTICAL_PARENT_SEQ,          
			 
			 
			 "h" => \$help,
	);

if ($help) { die $usage; }

unless ($query_NAST || $query_FASTA) { die $usage; }



main: {
	
	if ($exec_dir) {
		chdir $exec_dir or die "Error, cannot chdir to $exec_dir";
	}
	
	
	## ensure that ncbi blast formatted database exists.
	unless (-e "$db_FASTA.nhr" && -e "$db_FASTA.nin" && -e "$db_FASTA.nsq") {
		my $cmd = "formatdb -i $db_FASTA -p F 2>/dev/null";
		&process_cmd($cmd);
	}


	if ($query_FASTA && ! $query_NAST) {

		$query_NAST = $query_FASTA . ".NAST";
		
		my $cmd = "$FindBin::Bin/../NAST-iEr/run_NAST-iEr.pl --query_FASTA $query_FASTA "
			. " --db_NAST $db_NAST "
			. " --db_FASTA $db_FASTA "
			. " > $query_NAST ";
		
		&process_cmd($cmd);
	}
		

	## Run ChimeraParentSelector
	
	my $baseOut = basename($query_NAST);

	my $CPS_outfile = $baseOut . ".CPS";
	
	my $cmd = "$FindBin::Bin/ChimeraParentSelector/run_chimeraParentSelector.pl "
		. " --query_NAST $query_NAST "
		. " --db_NAST $db_NAST "
		. " --db_FASTA $db_FASTA ";

	$cmd .= " -n $numSeqsCompare " if defined($numSeqsCompare);
	$cmd .= " -P $minPerID " if defined($minPerID);
	$cmd .= " -R $minDivR " if defined($minDivR);
	
	$cmd .= " -M $matchScore " if defined($matchScore);
	$cmd .= " -N $mismatchPenalty " if defined($mismatchPenalty);
	$cmd .= " -T $maxTraverse " if defined($maxTraverse);
	$cmd .= " -Q $minQueryCoverage " if defined($minQueryCoverage);

	$cmd .= " --ignore_parent_same_acc_substring " if ($IGNORE_PARENT_SAME_ACC_SUBSTRING);
	$cmd .= " --ignore_identical_parent_seq " if ($IGNORE_IDENTICAL_PARENT_SEQ);
	
	$cmd .= " > $CPS_outfile";
	
	
	&process_cmd($cmd); 
	

	#############################
	####  RE-NAST align based on CPS results

	my $query_RENAST_file = "$query_NAST.CPS_RENAST";
	
	$cmd = "$FindBin::Bin/ChimeraParentSelector/CPS_to_RENAST.pl "
		. "--CPS_output $CPS_outfile "
		. "--query_NAST $query_NAST "
		. "--db_NAST $db_NAST "
		. "> $query_RENAST_file";
	
	&process_cmd($cmd); 
	
	


	########################
	#### Run ChimeraPhyloChecker
	
	my $CPC_outfile = $CPS_outfile . ".CPC";
		
	## Run ChimeraSlayer
	$cmd = "$FindBin::Bin/ChimeraPhyloChecker/CPS_to_CPC.pl "
		. " --CPS_output $CPS_outfile "
		. " --query_NAST $query_RENAST_file "
		. " --db_NAST $db_NAST ";
	
	
	$cmd .= " --num_parents_test $num_parents_test " if defined($num_parents_test);
	$cmd .= " -P $minPerID " if defined($minPerID);
	$cmd .= " -R $minDivR " if defined($minDivR);
	$cmd .= " -S $minPercentSampleSNPs " if defined($minPercentSampleSNPs);
	$cmd .= " --minBS $minBS " if defined ($minBS);
	
	$cmd .= " --num_BS_replicates $num_BS_replicates " if defined($num_BS_replicates);

	$cmd .= " --low_range_finer_BS $low_range_finer_BS " if defined($low_range_finer_BS);
	$cmd .= " --num_finer_BS_replicates $num_finer_BS_replicates " if defined($num_finer_BS_replicates);
	
	
	$cmd .= " --windowSize $windowSize " if defined($windowSize);
	$cmd .= " --windowStep $windowStep " if defined($windowStep);

	$cmd .= " --printAlignments " if ($printCSalignments);
	
	$cmd .= " > $CPC_outfile ";
	
	&process_cmd($cmd);
	
	## Add species and chimera type information
	$cmd = "$FindBin::Bin/util/CS_add_taxonomy.pl < $CPC_outfile > $CPC_outfile.wTaxons";
	&process_cmd($cmd);
	

	if ($printFinalAlignments) {
		## Run ChimeraParentSelector again to print alignments and reestimate breakpoints
		my $outfile = "$CPC_outfile.align";
		
		$cmd = "$FindBin::Bin/ChimeraParentSelector/CPC_to_CPS.pl "
			. " --CPC_output $CPC_outfile "
			. " --query_NAST $query_RENAST_file "
			. " --db_NAST $db_NAST "
			. " > $outfile ";
		
		&process_cmd($cmd);
	}
	
	exit(0);
}


####
sub process_cmd {
	my ($cmd) = @_;
	
	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd ($cmd) died with ret($ret)";
	}

	return;
}
