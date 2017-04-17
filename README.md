# rf_distances_filter
Scripts for analyzing Robinson-Foulds distances from Simmons et al. 2016

September 11, 2015
Dan Sloan
dbsloan@rams.colostate.edu

This directory contains three scripts used to compare and organize sets of gene trees based on Robinson-Foulds (RF) distances. There are two Python scripts (RFdistances.py and RFdistances.twofiles.py) that are modified versions of a script written by Sabrina Pankey (https://scriptomika.wordpress.com/2014/01/27/59/). RFdistances.py takes a set of tree files in newick format (one tree per file) as input and reports the RF distance and number of shared taxa for each pairwise comparison between trees. RFdistances.twofiles.py is similar except that the script only takes two tree files as input but each file can contain multiple trees. It returns the RF distance and number of shared taxa for each possible pairwise comparison involving a tree in file 1 and a tree in file 2. There is also a Perl wrapper script (RFdistances.filter.pl) that is used to call RFdistances.py and to return subsets of trees based on the resulting output from RFdistances.py. More details on each script are provided below. The directory also contains a folder of sample input tree files.


General Information and Dependencies

These scripts have been successfully installed on machines running Mac OSX 10.10 with Python 2.7 and Perl 5.18. They have not been tested for broader compatibility.

The python scripts require the following modules:

dendropy (version 3): https://github.com/jeetsukumaran/DendroPy/releases/tag/v3.12.1
p4 :  http://p4.nhm.ac.uk/



RFdistances.py

The RFdistances.py script would typically be called from the command line as follows:

python RFdistances.py -o myOutputFile.txt -p treeFilePrefix

In this case, myOutputFile.txt could be any output file name you choose and treeFilePrefix is a string found at the beginning of all of your desired input tree files (e.g. RAxML_result). The -p argument should include relevant path information if the tree files are not in the present working directory.

The output file contains four tab-delimited columns: 1) first tree; 2) second tree; 3) RF distance between the two trees; 4) number of shared taxa between the two trees.


RFdistances.twoFiles.v2.py

The RFdistances.twoFiles.v2.py script would typically be called from the command line as follows:

python RFdistances.twoFiles.v2.py -a File1 -b File2 -o myOutputFile.txt

File1 and File2 refer to the two tree files (newick format), which can (but do not have to) contain multiple trees each. The output file name can be any that you choose. The output format is the same as for RFdistances.py.



RFdistances.filter.pl


The usage statement for RFdistances.filter.pl is as follows:


Usage: RFdistances.filter.pl [options] treeFilePrefix 

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
RFdistances.filter.pl --subsets=80,85,90,95 --output=myOutput RAXML_bestTree


This script calls RFdistances.py, and the treeFilePrefix is the option that will be passed to RFdistances.py with the -p option.

After running RFdistances.py, this script will analyze the output and sort trees based on their average RF distances with the rest of the trees in the dataset. It returns subsets of trees with the lowest average RF distances. By default, it will return subsets with the lowest 10%, 20%, 30%, etc. But the cutoffs for these subsets can also be specified manually with the --subsets option. If the --rand_sets options is specified, the script will also generate random samples of trees with the same set sizes as defined by the --subsets option. 






