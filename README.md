# Bacillus project
Genome assembly, annotation and comparative genomics of 25 Bacillus species.

## Scope
zitw o trump!
## Assembly
The new species were assembled using the Perfect Genome Assembly tutorial.

How many isolates?

What sequencing technologies?

### Assembly pipeline

## Phylogeny
We used the `gtdb-tk` we build the phylogenetic tree of the assemblies.

First we run the `classify_wf` command to examine which assemblies can be 
classified based on the ANI. There were 3 assemblies that didn't classify
and some that are classified to the same species but we know from the lab
that they exhibit different properties.

The three unclassified assemblies we further explored with the `de_novo_wf` command,
in order to infer their phylogenetic similarity with other species in the tree.

The taxonomy was transformed with the following oneliner:

```
gawk -F"\t" 'BEGIN{print "gtdb_id" "\t" "classification_full" "\t" "classification" "\t" "taxonname"}{taxon=$1; split($2,taxonomy, "; ") ; for (i in taxonomy){split(taxonomy[i], classification,"__"); print taxon "\t" taxonomy[i] "\t" classification[1] "\t" classification[2]}}' gtdbtk.bac120.decorated.tree-taxonomy > gtdbtk.taxonomy_long.tsv
```
Quick look whether the assembly is in the tree.
```
gawk '{match($0,/543_assembly/,arr)}END{for (i in arr){print arr[i] "\t" arr[i,"start"] "\t" arr[i, "length"]}}' gtdbtk.bac120.decorated.tree
```
Then the results are loaded in R.

## Coding environment

Most computations were performed in the SarrisLab server: Linux sarrislab 22.04.1-Ubuntu,
Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz, 24 CPU cores, 256 gb Memory.

Additional computations were performed in the Zorbas HPC facility of [IMBBC-HCMR](https://hpc.hcmr.gr),
see here for more [info](https://doi.org/10.1093/gigascience/giab053).


