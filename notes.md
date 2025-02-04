# Genome annotation using Prokka

## Prokka Installation

I will install the Prokka package through bioconda using the following command:

```
conda create --name prokka
```
```
conda activate prokka
```
```
conda install bioconda::prokka
```

Update prokka:

```
conda install -c conda-forge -c bioconda prokka=1.14.5
```

## Use Prokka for the annotation of the assemblies

```
prokka --outdir /home/nik_arapitsas/Documents/Bacillus_project/Results/prokka_annotation --prefix SRL152_prokka_annotation /mnt/assemblies_repository/SRL152_assembly/SRL152_assembly.fasta 
```

```
art /home/nik_arapitsas/Documents/Bacillus_project/Results/prokka_annotation/SRL152_prokka_annotation.gff
```

# Comparative genomics (All vs All) using Orthofinder

## Unzip the IMG protein files from the zip files

```
mkdir -p /mnt/assemblies_repository/proteins_Bacillus_project  
```

The extract_proteins.sh script was created for the extraction of the protein files of the 25 Bacillus, from the zip files downloaded from the IMG database.

```
chmod +x extract_proteins.sh
```

```
./extract_proteins.sh
```

## Install Orthofinder

1) **Create an environment for Orthofinder**

```
conda create --name orthofinder
conda activate orthofinder
```

2) **Install Orthofinder**

```
conda install bioconda::orthofinder
```

3) **Update Orthofinder**

```
conda install -c conda-forge -c bioconda orthofinder=2.5.5
```

## Run Orthofinder

Diamond version was outdated and that made an issue in running Orthofinder. Thus it got updated to the last version sing the following command: 

```
conda update -c bioconda diamond
```

```
orthofinder -f /mnt/assemblies_repository/proteins_Bacillus_project -t 20 -o /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder
```

It took about 15 minutes to complete. It produced 741.5 Mb of output files. 

## Create Graphs

```
mkdir /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs
```

1) **Number of Isolate-Specific Orthogroups per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/species_specific_orthogroups.txt  
```

The graph was designed in Rstudio. The script is located in the following path:

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/orthofinder_graphs.R
```

2) **Percentage of genes from each isolate assigned to orthogroups**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

The same sorting with the list above will be used to provide plots that could be easily comparable:

First, extract the species order from the previous output file with species-specific orthogroups:

```
awk '{print $1}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/species_specific_orthogroups.txt > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/isolates.txt 
```

Then use it to assign the percentage of genes in orthofroups in the desirable order: 

```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/percofgenes_inogs_per_isolate_unsorted.txt /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/isolates.txt > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/percofgenes_inogs_per_isolate.txt
```

