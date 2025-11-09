# "Unravelling the genomic and functional arsenal of Bacilli endophytes from plants with different lifestyles"

This repository documents the computational pipeline used to assemble, annotate, and do comparative genomic analyses on the genomes of 25 Bacilli endophytes.

It aims to enhance the transparency and reproducibility of the pipeline used in our work with the title "Unravelling the genomic and functional arsenal of Bacilli endophytes from plants with different lifestyles".

## Genome assembly and annotation

### Filtering of Illumina short reads

Illumina short reads were filtered with [fastp v0.23.4](https://github.com/OpenGene/fastp).

The fastp filtering of the isolate SRL662 is provided below as an example:

```
fastp -i A01_FDSW210370227-1r_HLG2FDSX2_L1_1.fq.gz -I A01_FDSW210370227-1r_HLG2FDSX2_L1_2.fq.gz -o A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -O A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --report_title "SRL662 fastp report" --unpaired1 A01_FDSW210370227-1r_HLG2FDSX2_L1_1_unpaired.fq.gz --unpaired2 A01_FDSW210370227-1r_HLG2FDSX2_L1_2_unpaired.fq.gz
```

### Filtering of PacBio long reads

The raw sequencing data were obtained in bam files that were converted to fq files using [samtools](https://github.com/samtools/samtools) as shown below: 

```
samtools fastq output.bc1011_1--bc1011_1.subreads.bam > A01_long_unfiltered.fastq
```

Then the PacBio long reads were filtered with [Filtlong v0.2.1](https://github.com/rrwick/Filtlong) using settings “--min_length 1000” and “--keep_percent 95”:

```
filtlong --min_length 1000 --keep_percent 95 A01_long_unfiltered.fastq > A01_long_filtered.fastq
```

### Hybrid assembly

[Unicycler v0.5.0](https://github.com/rrwick/Unicycler) with its default settings was used for the hybrid assemblies. For the isolates **SRL389**, **SRL543** and **SRL662**, a first assembly was produced with [flye v.2.9.2-b1786](https://github.com/mikolmogorov/Flye). Then, this assembly was used in Unicycler with the --existing_long_read_assembly option, using the filtered Illumina and PacBio reads.

Below we will show the procedure combining flye and unicycler for the hybrid genome assembly of the isolate **SRL662** as an example. For the isolates where unicycler was performed directly, without a prior flye assembly, the --existing_long_read_assembly option in the unicycler command was omitted.    

#### A first long-read-only assembly using flye 

The command used for the isolate **SRL662** is the following:

```
flye --pacbio-raw A01_long_filtered.fastq --out-dir ../SRL662_flye_assembly --threads 23
```

#### Hybrid assembly using Unicycler with the existing long-read assembly parameter  

```
unicycler -1 A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l A01_long_filtered.fastq --existing_long_read_assembly ../SRL662_flye_assembly/assembly.fasta -o ../SRL662_flye_hybrid_assembly --threads 23
```

### QUAST and BUSCO analysis

The structural statistics of the assemblies were evaluated using [QUAST v5.0.2](https://github.com/ablab/quast) with default options: 

```
quast assembly.fasta -o ../SRL662_quast
```

For evaluating functional completeness, [BUSCO 5.8.0](https://github.com/metashot/busco) was used, with the Firmicutes database (firmicutes_odb10):

```
busco -i assembly.fasta -o SRL662_busco --out_path ../SRL662_busco -l firmicutes_odb10 -m genome -c 23 -f
```

### Annotation with Bakta

The annotation of every genome was generated using [Bakta v1.11.2](https://github.com/oschwengers/bakta) with its [full database v6.0](https://doi.org/10.5281/zenodo.14916843). 

#### Bakta installation using conda

```
conda create -n bakta
conda install -c conda-forge -c bioconda bakta
```

The full Bakta database was downloaded using the following command:

```
bakta_db download --output  --type full
```

** The databse had an unziped size of 89.6 Gb

#### Bakta annotation

The bakta command for SRL662 is presented below:

```
bakta --db db --prefix SRL662_bakta --output ../SRL662_bakta --threads 23 assembly.fasta
```

The [run_bakta.sh](scripts/run_bakta.sh) script was used for automatically performing bakta annotation to every assembly. 

### Plasmid prediction with RFPlasmid

[RFPlasmid v0.0.18](https://github.com/aldertzomer/RFPlasmid) was used to provide indications for the presence of plasmids in our assemblies using the [run_rfplasmid.sh](scripts/run_rfplasmid.sh) script.

## Taxonomic and phylogenetic analyses

The pipeline followed for the taxonomic and phylogenetic analyses was developed by Mrs Franziska Reden and Prof. Alexandros Stamatakis, and can be found in a separate [repository](https://github.com/FranziskaReden/bacillus_project). 

## OrthoFinder analysis 

The predicted proteomes of the isolates obtained from the Bakta annotation were collected in the same directory with the [collect_bakta_proteins.sh](scripts/collect_bakta_proteins.sh) script.  

Then, [OrthoFinder v.2.5.5](https://github.com/davidemms/OrthoFinder) was performed with default parameters as shown below:

```
orthofinder -f ./Bacillus_project_proteins -t 23 -o ./Bacillus_project_orthofinder
```

### Collection of the Data plotted on the collective OrthoFinder plot  

The OrthoFinder produced the drectory: 

```
Bacillus_project/Bacillus_project_orthofinder/Results_Jul03
```

The name of this directory changes in regard to the date when the OrthFinder command is performed by the user. We created a directory named "Graphs" in the directory "Results_Jul03" where we saved the data we used for the OrthoFinder graphs and all the outputs directly and indirectly associated with the graphs. 

```
cd Bacillus_project_orthofinder/Results_Jul03
mkdir Graphs
```

#### 1) Number of Isolate-Specific Orthogroups per Isolate

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > Graphs/species_specific_orthogroups.txt 
```

#### 2) Percentage of genes from each isolate assigned to orthogroups

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

In order to have in the list of the percentage of genes assigned to orthogroups for each isolates, the same sorting with the list of the number of species-specific orthogroups, we firstly extracted the isolates order from the previous output file (species_specific_orthogroups.txt):

```
awk '{print $1}' Graphs/species_specific_orthogroups.txt > Graphs/isolates.txt 
```

Then we used it to assign the percentage of genes in orthofroups in the same order: 

```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' Graphs/percofgenes_inogs_per_isolate_unsorted.txt Graphs/isolates.txt > Graphs/percofgenes_inogs_per_isolate.txt
```

#### 3) Orthogroups shared across all or any isolates

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    species[i] = substr($i, 1, 6); # Keep only first 6 letters
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  if (core) {
    for (i=2; i<=NF-1; i++) core_count[i]++; # Count core orthogroups for each species
    next; # Skip counting this orthogroup in the shared category
  }

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break}
      }
    }
  }
}
END {
  print "Isolates" "\t" "Core Orthogroups" "\t" "Partially Shared Orthogroups";
  for (i in species) {
    print species[i] "\t" core_count[i] "\t" shared[i];
  }
}
' Orthogroups/Orthogroups.GeneCount.tsv > Graphs/orthogroupcount_in_isolates.txt
```

The collective OrthoFinder graph was designed using the [orthofinder_graphs.R](scripts/orthofinder_graphs.R) script. 

## Pangenome analysis with anvi'o

[Anvi'o v8](https://github.com/merenlab/anvio) was used for the pangenome analysis. 

### Get the publicly available genomes of the close relatives of each isolate 

The publicly available genomes of the close relatives of each isolate was downloaded using the [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) command (v0.3.3). 

#### Make the directories where the genomes will be downloaded 

```
mkdir Bacillus_project_anvio
mkdir Bacillus_project_anvio/SRL179_anvio
mkdir Bacillus_project_anvio/SRL179_anvio/SRL179_genomes
mkdir Bacillus_project_anvio/SRL337_anvio
mkdir Bacillus_project_anvio/SRL337_anvio/SRL337_genomes
mkdir Bacillus_project_anvio/SRL368_anvio
mkdir Bacillus_project_anvio/SRL368_anvio/SRL368_genomes
mkdir Bacillus_project_anvio/SRL543_anvio
mkdir Bacillus_project_anvio/SRL543_anvio/SRL543_genomes
```

#### i) SRL179

```
ncbi-genome-download bacteria --assembly-accessions GCF_030123405.1,GCF_000612625.1,GCF_023715505.1 --formats fasta --flat-output --output-folder ./Bacillus_project_anvio/SRL179_anvio/SRL179_genomes 
gunzip ./Bacillus_project_anvio/SRL179_anvio/SRL179_genomes/*.fna.gz
```

#### ii) SRL337

```
ncbi-genome-download bacteria --assembly-accessions GCF_003581585.1,GCF_008180615.1,GCF_008180415.1 --formats fasta --flat-output --output-folder ./Bacillus_project_anvio/SRL337_anvio/SRL337_genomes 
gunzip ./Bacillus_project_anvio/SRL337_anvio/SRL337_genomes/*.fna.gz
```

#### iii) SRL368

```
ncbi-genome-download bacteria --assembly-accessions GCF_020775355.1,GCF_022663675.1,GCF_000021305.1,GCF_002559825.1,GCF_001182785.1 --formats fasta --flat-output --output-folder ./Bacillus_project_anvio/SRL368_anvio/SRL368_genomes 
gunzip ./Bacillus_project_anvio/SRL368_anvio/SRL368_genomes/*.fna.gz
```

#### iv) SRL543

```
ncbi-genome-download bacteria --assembly-accessions GCF_024706405.1,GCF_020171945.1,GCF_000473245.1 --formats fasta --flat-output --output-folder ./Bacillus_project_anvio/SRL543_anvio/SRL543_genomes 
gunzip ./Bacillus_project_anvio/SRL543_anvio/SRL543_genomes/*.fna.gz
```

**There was a newer version for the assembly GCF_024706405.1, the GCF_024706405.2**. So we downloaded this:

```
ncbi-genome-download bacteria --assembly-accessions GCF_024706405.2 --formats fasta --flat-output --output-folder ./Bacillus_project_anvio/SRL543_anvio/SRL543_genomes 
gunzip ./Bacillus_project_anvio/SRL543_anvio/SRL543_genomes/*.fna.gz
```

### Pangenome analysis 

The pipeline for SRL179 pangenome analysis with its closest relatives will be presented as an example. The same procedure was followed for the rest pangenome analyses. 


1) Rename the files, change the contig naming in the fasta files and prepare a database for every genome

```
cd ./Bacillus_project_anvio/SRL179_anvio
mkdir ./SRL179_genomes_simplified 
mkdir ./SRL179_genomes_db
cd ./SRL179_genomes
```

```
for file in *.fna *.fasta; do
  [ -e "$file" ] || continue
  if [[ "$file" == GCF_*_genomic.fna ]]; then
    # Extract GCF prefix and version, replace '.' with '_'
    newname=$(echo "$file" | sed -E 's/^(GCF_[0-9]+)\.([0-9]+).*\.fna$/\1_\2.fna/')
    mv "$file" "$newname"
  elif [[ "$file" == SRL* ]]; then
    newname=$(echo "$file" | sed -E 's/^(SRL[0-9]+).*\.fasta$/\1.fasta/')
    mv "$file" "$newname"
  fi
done
```

```
for fasta in *.fasta *.fna; do
  [ -e "$fasta" ] || continue
  base=$(basename "$fasta")
  prefix="${base%.*}"

  # Reformat fasta and simplify contig names
  anvi-script-reformat-fasta "$fasta" \
    -o "../SRL179_genomes_simplified/${prefix}_simplified.fna" \
    --seq-type NT \
    --simplify-names \
    --report-file "../SRL179_genomes_simplified/${prefix}_rename-report.txt" \
    --prefix "$prefix"
  
  anvi-gen-contigs-database -f "../SRL179_genomes_simplified/${prefix}_simplified.fna" -o "../SRL179_genomes_db/${prefix}.db" --project-name "$prefix" -T 20 

done
```

2) Annotate the db files

```
anvi-setup-kegg-data
anvi-setup-ncbi-cogs
anvi-setup-pfams
anvi-setup-scg-taxonomy
anvi-setup-trna-taxonomy
anvi-setup-interacdome
```
```
cd ../SRL179_genomes_db 
```
```
for db in *.db; do
  echo "🔍 Annotating $db"
  echo "🔍 1. HMMS Annotation"
  anvi-run-hmms -c "$db" -I Bacteria_71 --also-scan-trnas -T 20
  echo "🔍 2. NCBI COG Annotation"
  anvi-run-ncbi-cogs -c "$db" -T 20
  echo "🔍 3. KEGG Kofam Annotation"
  anvi-run-kegg-kofams -c "$db" -T 20
  echo "🔍 4. Pfam Annotation"
  anvi-run-pfams -c "$db" -T 20
  echo "🔍 5. SCG Annotation"
  anvi-run-scg-taxonomy --contigs-db "$db" -T 20
done
```

3) Create a genomes-storage database file

```
echo -e "name\tcontigs_db_path" > external_genomes.txt
for db in *.db; do
  name=$(basename "$db" .db)
  echo -e "${name}\t$(realpath "$db")" >> external_genomes.txt
done
```
```
anvi-gen-genomes-storage -e external_genomes.txt -o SRL179-GENOMES.db
```

4) Run pangenome analysis

```
anvi-pan-genome -g SRL179-GENOMES.db \
                --project-name "SRL179_Pangenome" \
                --output-dir ../SRL179_pangenome \
                --num-threads 20 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
```

### Displaying the pan genome

```
anvi-display-pan -p ../SRL179_pangenome/SRL179_Pangenome-PAN.db -g ./SRL179-GENOMES.db
```

### Extraction of gene cluster presence/absence data from the pangenome analysis

```
cd ../SRL179_pangenome
```

```
sqlite3 -header -separator $'\t' SRL179_Pangenome-PAN.db "SELECT gene_cluster_id, genome_name FROM gene_clusters;" > SRL179_gene_cluster_output.txt
```
```
gawk -F"\t" '(NR>1){a[$1FS$2]++}END{for (i in a){print a[i] FS i}}' SRL179_gene_cluster_output.txt | gawk -F"\t" '{a[$2]++}END{for (i in a){print a[i] FS i}}' | sort -k1 -r | head 
```

### Design of the Pangenome Graphs

A directory for the pangenome graphs was created in the "Bacillus_project_anvio" directory, using the following command:

```
mkdir Bacillus_project_anvio_Graphs
```

The files SRL179_gene_cluster_output.txt, SRL337_gene_cluster_output.txt, SRL368_gene_cluster_output.txt and SRL543_gene_cluster_output.txt were collected inside the "Bacillus_project_anvio_Graphs". Then, the graphs were created in R (version 4.5.1) using the [pangenome_graphs.R](scripts/pangenome_graphs.R) script. This script was executed inside the "Bacillus_project_anvio_Graphs" directory.  


### Biosynthetic Gene Clusters (BGCs) analysis

For BGC prediction, the assemblies were uploaded to the online platform of [antiSMASH version 8.0.1](https://antismash.secondarymetabolites.org/#!/start).   

The information regarding the BGCs of each isolate and their similarity confidence with known BGCs was manually collected from the output of antismash. This information was used to prepare the tables [antismash_table.csv](data/antismash_table.csv) and [bgc_similarity_with_species.csv](data/bgc_similarity_with_species.csv) that are available in the [data directory](data) of this repository. 

Graphs were prepared in R (version 4.5.1) using the [antismash_bgcs_pie_charts.R](scripts/antismash_bgcs_pie_charts.R), [antismash_groups_similarity_graphs.R](scripts/antismash_groups_similarity_graphs.R), and [bgcs_barplot_pergenus.R](scripts/bgcs_barplot_pergenus.R) scripts.  


### Plotting of gene counts associated with plant–microbe interactions across the 25 bacterial isolates using a heatmap

Gene counts across nine functional categories and isolates were collected in the file [Functions_Distribution_per_Isolate.xlsx](data/Functions_Distribution_per_Isolate.xlsx) that is available in the [data directory](data) of this repository. For their visualization, the script [gene_function_heatmap.R](scripts/gene_function_heatmap.R) was used in R (version 4.5.1). 


## Coding environment

Most computations were performed in the SarrisLab server: Linux sarrislab 22.04.1-Ubuntu,
Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz, 24 CPU cores, 256 gb Memory.

Additional computations were performed in the Zorbas HPC facility of [IMBBC-HCMR](https://hpc.hcmr.gr),
see here for more [info](https://doi.org/10.1093/gigascience/giab053).



