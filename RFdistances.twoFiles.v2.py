#! /usr/bin/env python
import p4, dendropy
import sys, glob
import types
import ast
import exceptions
from optparse import OptionParser
import cStringIO
import copy
 
def setuploop(file_1,file_2,outfile):
	p4.var.trees = []
	p4.read(file_1)
	fileCount1 = len(p4.var.trees)
	p4.read(file_2)
	fileCount2 = len(p4.var.trees) - fileCount1
	
	for i in range(0, fileCount1):
		for j in range(fileCount1, fileCount1+fileCount2):
			
			temp_tree1 = copy.deepcopy(p4.var.trees[i])
			temp_tree2 = copy.deepcopy(p4.var.trees[j])
			
			tree1,tree2,numberOfTaxa = _pruneTreesToComplement(temp_tree1, temp_tree2)
			if tree1 != 'nope':
				t1 = dendropy.Tree.get_from_string(tree1, 'newick')
				t2 = dendropy.Tree.get_from_string(tree2, 'newick')
				symdiff = t1.symmetric_difference(t2) #/ len(t1.nodes())
				towrite2 = file_1 + '\t' + file_2 +'\t' + str(symdiff) + '\t' + str(numberOfTaxa) + '\n'
				with open(outfile, 'a') as out2:
					out2.write(towrite2)
         
def _pruneTreesToComplement(tree1,tree2):
 
   # Delete out the missing taxa from trees so that they intersect:
    tree1_taxa = [n.name for n in tree1.nodes if n.isLeaf]
    tree2_taxa = [n.name for n in tree2.nodes if n.isLeaf]
    commontaxa = list(set(tree1_taxa) & set(tree2_taxa))
    numberOfSharedTaxa = len(commontaxa)
    if len(commontaxa) > 3:
      delete_from_tree1 = set([t for t in tree1_taxa if t not in tree2_taxa])
      delete_from_tree2 = set([t for t in tree2_taxa if t not in tree1_taxa])
      for t in delete_from_tree1:
          tree1.removeNode(t, alsoRemoveBiRoot=False)
      for t in delete_from_tree2:
          tree2.removeNode(t, alsoRemoveBiRoot=False)
      if len([n.name for n in tree1.nodes if n.isLeaf]) != \
            len([n.name for n in tree2.nodes if n.isLeaf]):
          raise TCTError, 'Something unexpected went wrong with the taxon ' + \
                'pruning. After pruning trees have different number of taxa.'
      t1=tree1.writeNewick(toString=True).split('\n')[0]
      t2=tree2.writeNewick(toString=True).split('\n')[0]
      return (t1,t2,numberOfSharedTaxa)
    else:
      return('nope','nope','nope')
def main():
    msg   = 'usage: ./compare_distance.py -a <file 1> -b <file 2> -o <output>' + \
             '\n\nReads two trees (NEXUS and/or PHYLIP format),\n' +\
             'reciprocally prunes each tree of missing taxa ' + \
             'and returns the distance between each pair of trees' + \
             '\n\nDepends on Python2.6, dendopy and p4 modules.\nMust run from within directory containing treefiles.'
    usage = msg
    parser = OptionParser(usage)
    parser.add_option('-a', '--file1',dest='first_file', help='Name of treefile 1')
    parser.add_option('-b', '--file2',dest='second_file', help='Name of treefile 2')
    parser.add_option('-o', '--out',dest='output', help='Name for output containing distances between trees')
    (options,args) = parser.parse_args()
    setuploop(options.first_file,options.second_file,options.output)
     
 
if __name__ == "__main__":
    main()