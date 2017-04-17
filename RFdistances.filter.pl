#!/usr/bin/perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use List::Util qw(sum shuffle);
use Getopt::Long;
use File::Which;


#parse/define arguments
our $SUBSETS = "10,20,30,40,50,60,70,80,90,100";
our $OUTPUT;
our $PYTHON = which('RFdistances.py');
$PYTHON or $PYTHON = 'RFdistances.py';
our $RAND_SETS = 0;
our $RF_FILE;
parseArgs();
my $prefix = shift;
$OUTPUT or $OUTPUT = $prefix;

#define subsets based on comma-delintied string of percentages
my %subsetContents;
my @splitSubsets = split (/\,/, $SUBSETS);
foreach (@splitSubsets){
	unless (looks_like_number($_) && $_ > 0 && $_ <= 100){
		die ("\nERROR: Specified list of subsets should be comma-delimited and contain non-zero values less than or equal to 100. Offending value: $_\n\n");
	}
	$subsetContents{$_} = "";
}

#check that python script can be found
-e $PYTHON or die ("\nERROR: Could not find the python script $PYTHON\. The full path can be specified with the --pyth_script option.\n\n");

#check that RAND_SETS is a number
looks_like_number($RAND_SETS) or die ("\nERROR: If specified --rand_sets must be followed by a numerical value\n\n");

#call python script unless pairwise RF file has been provided as input
unless ($RF_FILE){ 
	system ("python $PYTHON -p $prefix -o $OUTPUT\.pairwiseRF.txt");
	$RF_FILE = "$OUTPUT\.pairwiseRF.txt";
}

#read in output from python script
my @lines = file_to_array ("$RF_FILE");

#define data structures to store calculated distances and filename info
my %RF_dist_HOA;
my %RF_avgDist;
my @sortedFiles;

#loop over each line in the python script output file
foreach (@lines){
	$_ =~ /^\s*$/ and next;
	chomp $_;
	
	my @sl = split (/\t/, $_);
	
	#calculate corrected RF distance as the raw distance divided by 2 * (number of taxa - 3)
	my $RFcorr = $sl[2]/(2*($sl[3] - 3));
	
	#add corrected RF distance to an array of distances for each tree1 in the input file
	#There will be a separate array for each different tree in column 1, and these are organized in a hash of arrays (tree file names as keys, and refs to distance arrays as values)
	push (@{$RF_dist_HOA{$sl[0]}}, $RFcorr);
}

#Calculate average corrected RF distances for each tree
foreach (keys %RF_dist_HOA){
	my @distArray = @{$RF_dist_HOA{$_}};
	my $avgDist = sum(@distArray)/@distArray;
	$RF_avgDist{$_} = $avgDist;	
}

#sort by avg RF value and print tab-delimited output
my $FH1 = open_output("$OUTPUT\.avg_dist.txt");
foreach (sort { $RF_avgDist{$a} <=> $RF_avgDist{$b} or $a cmp $b } keys %RF_avgDist) {
	print $FH1 "$_\t$RF_avgDist{$_}\n";
	push (@sortedFiles, $_);
}
close $FH1;


#concatenate trees into different subsets based on sorted distances (e.g., the subset 90 includes the 90% of trees with the lowest average RF distances)
my $totalTrees = scalar (keys %RF_avgDist);
my $treeNum = 0;
foreach my $treeFile (@sortedFiles){
	++$treeNum;
	my $fileContents = file_to_string ($treeFile);
	foreach (keys %subsetContents){
		if ($treeNum/$totalTrees <= $_/100){
			$subsetContents{$_} .= $fileContents;
		}
	}	
}

#print concated tree output to files
foreach (keys %subsetContents){
	my $FH2 = open_output ("$OUTPUT\.subset$_\.tre");
	print $FH2 $subsetContents{$_};
}

#generate and print random subsets (as above but with a randomly shuffled list of files)
if ($RAND_SETS){
	my $rand_count = 1;
	while ($rand_count <= $RAND_SETS){
		my @randFiles = shuffle @sortedFiles;
		my $randTreeNum = 0;
		my %randSubsetContents;
		foreach my $randTreeFile (@randFiles){
			++$randTreeNum;
			my $fileContents = file_to_string ($randTreeFile);
			foreach (keys %subsetContents){
				if ($randTreeNum/$totalTrees <= $_/100){
					$randSubsetContents{$_} .= $fileContents;
				}
			}	
		}
		foreach (keys %randSubsetContents){
			my $FH3 = open_output ("$OUTPUT\.rand$rand_count\.subset$_\.tre");
			print $FH3 $randSubsetContents{$_};
		}			
		++$rand_count;
	}	
}


sub parseArgs{

  my $usage = 

"\nUsage: $0 [options] treeFilePrefix 
   
   ARGUMENTS
         treeFilePrefix - Text shared at the beginning of all input tree files
                          (Include path if necessary)
         
   	
   OPTIONS
         --subsets      - a comma-delimited list of % cutoffs for subet files
                          default: 10,20,30,40,50,60,70,80,90,100

         --output       - a base name for all output files generated
                          (Include path if necessary)
                          default: treeFilePrefix

         --pyth_script  - full path to RFdistances.py script
                          default: RFdistances.py
                          
         --rand_sets    - number of random subset files to generate                 
                          default: 0
                          
         --RF_file      - full path to file containing output from pairwise RF
                          calculations to avoid re-running this long step                 
                          default: 0


   EXAMPLE
         $0 --subsets=80,85,90,95 --output=myOutput RAXML_bestTree
\n\n";
  
  my $result = GetOptions
    (
     'subsets=s'  => \$SUBSETS,
     'output=s'  => \$OUTPUT,
     'pyth_script=s'  => \$PYTHON,
     'rand_sets=s'  => \$RAND_SETS,
     'RF_file=s'  => \$RF_FILE,
    );
  
  unless( @ARGV == 1 ){
    print $usage and exit;
  }
  
}



# A subroutine to get data from a file given its filename and store as array
sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

# A subroutine to get data from a file given its filename and store as string
sub file_to_string {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my $filedata;

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    
    while (<GET_FILE_DATA>){
    	$filedata .= $_;
    }
    
    close GET_FILE_DATA;

    return $filedata;
}

#A subroutine to open output file and return file handle
sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}
