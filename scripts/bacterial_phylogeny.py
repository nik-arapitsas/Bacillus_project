#!/usr/bin/env python3

###############################################################################
# script name: bacterial_phylogeny.py
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:

###############################################################################
# usage:./bacterial_phylogeny.py
###############################################################################
import os,sys
import dendropy

print("Loading the tree")

path = "/home1/s.paragkamian/bacillus/543_de_novo/"
filename = "gtdbtk.bac120.decorated.tree"
tax = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
    path=os.path.join(path, filename),
    schema="newick",
    suppress_internal_node_taxa=True,
    taxon_namespace=tax)

print(tree.seed_node)
print(type(tree.taxon_namespace))

for nd in tree:
    print(nd.label)
    print(nd.taxon)

print(len(tree.taxon_namespace))
taxa_to_retain = set([taxon for taxon in tree.taxon_namespace
        if taxon.label.startswith("f__DSM-18226")])
tree1 = tree.extract_tree_with_taxa(taxa=taxa_to_retain)

print(tree1.as_ascii_plot())

sys.exit()
print("Writing the tree")
output = "tree_ascii.txt"
tree_ascii = open(os.path.join(path, output),"w")

tree_ascii.write(tree.as_ascii_plot())
tree_ascii.close()


filter = lambda taxon: True if taxon=="543_assebly" else False
print(filter)
#node = tree.find_node_with_taxon("'f__DSM-18226'")
print(node)
# subtree 
#taxa_to_retain = {"f__DSM-18226"}
#tree2 = tree.extract_tree_with_taxa(taxa=taxa_to_retain)

#print(tree2.as_ascii_plot())
sys.exit()
#tree.taxon_namespace = taxa
#node = tree.ancestor_iter('RS_GCF_020171945.1')
#print(node)

#print("Print the distances")
#pdc = tree.phylogenetic_distance_matrix()
#for i, t1 in enumerate(tree.taxon_namespace[:-1]):
#    for t2 in tree.taxon_namespace[i+1:]:
#        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

#sys.exit()
