#! /usr/bin/env python
import p4, dendropy
import sys, glob
import types
import ast
import exceptions
from optparse import OptionParser
import cStringIO
 
def setuploop(treefilelist,outfile):
    for file in treefilelist:
        genename1 = file #saves last half of filename as gene. should match phytab
        position = treefilelist.index(file)
        allotherfiles = [f for f in treefilelist if treefilelist.index(f) != treefilelist.index(file)]
        for f in allotherfiles:
           genename_f = f
           tree1,tree2,numberOfTaxa = _pruneTreesToComplement(file, f)
 #use dendropy to get distance, save it, and then continue loop
           if tree1 != 'nope':
             t1 = dendropy.Tree.get_from_string(tree1, 'newick')             
             t2 = dendropy.Tree.get_from_string(tree2, 'newick') 
             symdiff = t1.symmetric_difference(t2) #/ len(t1.nodes())
             towrite2 = genename1 + '\t' + genename_f +'\t' + str(symdiff) + '\t' + str(numberOfTaxa) + '\n'
             with open(outfile, 'a') as out2:
                 out2.write(towrite2)
         
def _pruneTreesToComplement(file1,file2):
    p4.var.trees = []
 
    p4.read(file1)
    p4.read(file2)
    tree1 = p4.var.trees[0]
    tree2 = p4.var.trees[1]
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
    msg   = 'usage: ./compare_distance.py -p <file prefix> -o <output>' + \
             '\n\nReads two trees (NEXUS and/or PHYLIP format),\n' +\
             'reciprocally prunes each tree of missing taxa ' + \
             'and returns the distance between each pair of trees' + \
             '\n\nDepends on Python2.6, dendopy and p4 modules.\nMust run from within directory containing treefiles.'
    usage = msg
    parser = OptionParser(usage)
    parser.add_option('-p', '--trpref',dest='treeprefix', help='Common prefix shared by all tree files in directory. eg: RAxML_result')
    parser.add_option('-o', '--out',dest='output', help='Name for output containing distances between trees')
    (options,args) = parser.parse_args()
    treeprefix = options.treeprefix
    listoftrees = [i for i in glob.glob(treeprefix+"*")]
    setuploop(listoftrees,options.output)
     
 
if __name__ == "__main__":
    main()