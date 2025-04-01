#!/bin/bash

gtdbtk de_novo_wf \
	--batchfile data/bacillus_assembly_list.txt \
	--out_dir gtdbtk_de_novo_all \
	--bacteria \
	--cpus 20 \
	--outgroup_taxon p__Bacillota


