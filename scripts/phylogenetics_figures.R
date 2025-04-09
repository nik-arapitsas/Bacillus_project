#!/usr/bin/env Rscript

###############################################################################
# script name: phygenetics_figures.R
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:
# Aim of this script is to analyse and  visualise the de novo trees from the
# gtdb-tk pipeline.
###############################################################################
# usage:./phygenetics_figures.R
###############################################################################

library(tidyverse)
library(ape)
library(ggtree)
library(tidytree)

# read the taxonomy
# gawk -F"\t" '{taxon=$1; split($2,taxonomy, "; ") ; print "gtdb_id" "\t" "classification_full" "\t" "classification" "\t" "taxonname" ; for (i in taxonomy){split(taxonomy[i], classification,"__"); print taxon "\t" taxonomy[i] "\t" classification[1] "\t" classification[2]}}' gtdbtk.bac120.decorated.tree-taxonomy > gtdbtk.taxonomy_long.tsv

taxonomy <- read_delim( "../gtdbtk_de_novo_all/infer/gtdbtk.taxonomy_long.tsv",
                       delim="\t",
                       col_names=T)
# transform taxonomy to wide format
taxonomy_w <- taxonomy %>% group_by(gtdb_id, classification) %>%
    summarise(taxa=paste(taxonname, collapse = "; "),
              .groups = "drop") %>%  
    pivot_wider(names_from=classification, values_from=taxa)

# load tree
tree <- ape::read.tree(file = "../gtdbtk_de_novo_all/gtdbtk.bac120.decorated.tree")
tree_table <- read_delim( "../gtdbtk_de_novo_all/gtdbtk.bac120.decorated.tree-table", delim="\t")
tree_taxonomy <- as_tibble(tree) %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
    treeio::as.treedata()

# filter the tree

tree_bacilli <- as_tibble(tree) %>% filter(grepl("Bacilli",label))

# it doesn't have a node for Bacilli.
# most recent common ancestor (MRCA)
nodes_assemblies <- as_tibble(tree) %>%
    filter(grepl("^SRL",label)) #"337_assembly"

all_assemblies <- as_tibble(tree) %>%
    filter(grepl("^SRL",label)) #"337_assembly"
    #filter(label %in% c("179_assembly", "543_assembly","337_assembly")) #"337_assembly"

all_assemblies_m <- MRCA(tree,all_assemblies$node ) 
assemblies_mrca <- MRCA(tree,nodes_assemblies$node )

tree_assemblies <- extract.clade(tree, assemblies_mrca) %>%
    as_tibble() %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
    as.treedata()

all_assemblies <- extract.clade(tree, all_assemblies_m) %>%
    as_tibble() %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
    mutate(species=ifelse(s=="NA",label,s)) %>%
    as.treedata()

all_assemblies_taxonomy <- all_assemblies |>
    as_tibble() |>
    filter(grepl("^SRL",label))

#claude_mrca <- groupClade(tree_taxonomy, assemblies_mrca) %>%
#    filter(group==1) 

#write.tree(as.phylo(claude_mrca), file = "claude_mrca.tree", append = FALSE,
#           digits = 10, tree.names = FALSE)

# the node 82141 is the family
claude_f__DSM_18226 <- groupClade(as_tibble(tree), 82141) %>%
    filter(group==1) %>% 
    left_join(taxonomy_w,
              by=c("label"="gtdb_id")) %>% 
    as.treedata()

# visualisation

## default
tree_plot_mrca <- ggtree(tree_assemblies) + 
     geom_tree(aes(color=o, fill=o)) + 
     geom_tiplab(aes(as_ylab=TRUE, color=o), size=0.7) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = "bottom")
    
ggsave(plot=tree_plot_mrca,
       "../tree_plot_mrca_gtdb.pdf",
       device="pdf",
       height = 300,
       width=40,
       units="cm", limitsize = F)

tree_plot_mrcaC <- ggtree(tree_assemblies, layout="inward_circular") + 
     geom_tree(aes(color=o, fill=o)) + 
     geom_tiplab(aes(as_ylab=TRUE, color=o), size=0.5) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = "bottom")
    
ggsave(plot=tree_plot_mrcaC,
       "../tree_plot_mrca_gtdb_c.pdf",
       device="pdf",
       height = 300,
       width = 300,
       units="cm", limitsize = F)

tree_plot_mrca2 <- ggtree(tree_assemblies) + 
     geom_tiplab(size=2, aes(label=s, color=g)) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))

ggsave(plot=tree_plot_mrca2,
       "../tree_plot_mrca_s.pdf",
       device="pdf",
       height = 80,
       width=40,
       units="cm", limitsize = F)


tree_plot_mrca3 <- ggtree(all_assemblies) + 
     geom_tiplab(size=2, aes(label=species, color=g)) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))

ggsave(plot=tree_plot_mrca3,
       "../tree_plot_mrca_3.pdf",
       device="pdf",
       height = 100,
       width=40,
       units="cm", limitsize = F)

