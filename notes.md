# Organization of the isolates data

Initially, the assemblies with their quast output were saved in the directory "/media/sarlab/DATA/Bacillus_project" in folders named "*Isolate code*_assembly (were Isolate code: the SRL followed by the number of each isolate), and the Results in the directory "/home/nik_arapitsas/Documents/Bacillus_project/Results". The changes in the directory "/media/sarlab/DATA/Bacillus_project" were made mostly automatically as presented below, while most of the rearrangements from the "/home/nik_arapitsas/Documents/Bacillus_project/Results" directory were done manually (except of the antiSMASH output). The movement of the antiSMASH file and the protein files will be provided later.

## Rearrangements in the "/media/sarlab/DATA/Bacillus_project" directory

Go to the "/media/sarlab/DATA/Bacillus_project" directory:  

```
cd /media/sarlab/DATA/Bacillus_project
```

With the below code every folder will keep only the first six characters, which is the name of every isolate:

```
for d in */; do mv "$d" "${d:0:6}"; done
```

Then I used the code below to create a directory called "*Isolate code*_assembly" for every isolate and then move every assembly-associated file in it: 

```
for d in */; do 
  mkdir "$d"/"${d:0:6}"_assembly
  find "$d" -mindepth 1 -maxdepth 1 -type f -exec mv {} "$d"/"${d:0:6}"_assembly/ \;
done
```

* The rearrangement of the protein files is presented in last part of the section **"Comparative genomics (All vs All) using Orthofinder"**.

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

**Eventually, I did not use it.**

# Comparative genomics (All vs All) using Orthofinder

## Unzip the IMG protein files from the zip files

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

**In order for Orthofinder to run, the protein sequence files must be in te same directory. So the orthofinder had been run when the protein sequences were all in the directory "/mnt/assemblies_repository/proteins_Bacillus_project". 

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

3) **Genes with orthogroups in all or any isolates**

With the code below when counting the partially shared orthogroups we do not count the core orthogroups.

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
' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Orthogroups/Orthogroups.GeneCount.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/orthogroupcount_in_isolates.txt
```

With the code below when counting the partially shared orthogroups we count the core orthogroups as well. This is better for creating a bar plot where the bars of all and any orthogroups will be overlapping. 

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    species[i] = substr($i, 1, 6); # Keep only first 6 letters of species name
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      if (core) core_count[i]++;  # Count genes in core orthogroups
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break} # Count genes in shared orthogroups, including core
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
' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Orthogroups/Orthogroups.GeneCount.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/orthogroupcount_in_isolates.txt
```

# Comparative genomics (All vs All) using Orthofinder for the SRL368 and relatives

## Run Orthofinder

```
orthofinder -f /mnt/assemblies_repository/proteins_Bacillus_th_israel -t 20 -o /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives
```

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. & Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914

## Create Graphs

```
mkdir /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Graphs
```

1) **Number of Isolate-Specific Orthogroups per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, index($i, "_")-1)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Graphs/species_specific_orthogroups.txt  
```

2) **Number of Isolate-Specific Genes per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, index($i, "_")-1)} NR==10 {for(i=2; i<=NF; i++) print species[i], $i}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Graphs/species_specific_genes_number.txt  
```

3) **Percentage of genes from each isolate assigned to orthogroups**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, index($i, "_")-1)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

4) **Genes with orthogroups in all or any isolates**

With the code below when counting the partially shared orthogroups we count the core orthogroups as well. This is better for creating a bar plot where the bars of all and any orthogroups will be overlapping. 

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    for(i=2; i<=NF; i++) species[i]=substr($i, 1, index($i, "_")-1);
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      if (core) core_count[i]++;  # Count genes in core orthogroups
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break} # Count genes in shared orthogroups, including core
      }
    }
  }
}
END {
  print "Isolates" "\t" "Core Orthogroups" "\t" "Partially Shared Orthogroups";
  for (i=2; i<=NF; i++) {  # Iterate over the expected index range
    if (species[i] != "") {  # Ensure valid species name
      print species[i] "\t" (core_count[i] ? core_count[i] : 0) "\t" (shared[i] ? shared[i] : 0);
    }
  }
}
' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Orthogroups/Orthogroups.GeneCount.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Graphs/orthogroupcount_in_isolates.txt
```

## Check which genes are in the species specific orthogroups of SRL368:

```
mkdir /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes 
```
### 1) Get the species-specific Orthogroups for SRL368

```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignoring last column
{
  target = 5;  # Column for isolate E (adjust if needed)
  if ($target > 0) {  # Ensure isolate E has genes
    is_species_specific = 1;

    # Check if any other isolate (excluding column 6) has genes
    for (i=2; i<=NF-2; i++) {  # NF-2 to ignore the last column
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Orthogroups/Orthogroups.GeneCount.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_orthogroups.txt
```

Initially I did not get anything. I found out that it was searching in column 6 (total) as well. That was the reason for the issue, as the column 6 will have always a value > 0. So, I changed the code in order not to search in the column 6. 

### 2) Get the species-specific Genes for SRL368 by searching using species-specific Orthogroups

First I need to use tab as delimiter in the txt file: 

```
sed 's/ \+/	/g' species_specific_orthogroups.txt > species_specific_orthogroups_with_tabs.txt
```

Validate it with:

```
cat -T species_specific_orthogroups_with_tabs.txt
```

As ^I symbols appeared it validated that the tab was set as a delimiter successfully.

In order to get the isolate-specific genes I used the following command:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from species_specific_orthogroups.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $5;  # Extract the gene names from the fifth column
        print orthogroup, genes;  # Print orthogroup name and corresponding gene names
    }
}
' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_orthogroups_with_tabs.txt /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Orthogroups/Orthogroups.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_genes.txt
```

Get as output the genes in a list:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from species_specific_orthogroups_with_tabs.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $5;  # Extract the gene names from the fifth column
        split(genes, gene_array, ",");  # Split the gene names into an array
        for (i in gene_array) {
            print orthogroup, gene_array[i];  # Print orthogroup and gene in two columns
        }
    }
}
' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_orthogroups_with_tabs.txt /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/Orthogroups/Orthogroups.tsv > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_genes_list.txt
```

Extract the protein sequence of each isolate specific gene: 

```
./extract_isolate_specific_genes.sh
```

### 3) Run blastp for the sequences

First create a fasta file with the sequences:

```
awk '{print ">"$1"\n"$2}' /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/filtered_protein_sequences.txt > /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/filtered_protein_sequences.fasta
```

Run BLASTP against the nr (non-redundant protein sequences) database (the default of the online BLAST) for every sequence of our file: 

```
conda activate perfect_assembly
```

```
blastp -query /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/filtered_protein_sequences.fasta -db nr -out /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/blastp_results_filtered_protein_sequences.txt -evalue 1e-5 -num_threads 20 -outfmt 6
```
blastp -db nr -query proteins.fasta -remote -out /home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/result.txt 

## Movement of the protein files to meet the new server directory architecture

Move all the protein files in the respect directory:

```
cd /home/nik_arapitsas/Desktop/Test/proteins_Bacillus_project

for file in *; do
  prefix="${file:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}/${prefix}_proteins"
  if [ -d "$target_dir" ]; then
    mv "$file" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for file '$file'"
  fi
done
```
Move also the .zip files that had been downloaded from the IMG daabase and were saved in the /Downloads/Test directory.

```
cd /home/nik_arapitsas/Downloads/Test

for file in *; do
  prefix="${file:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}/${prefix}_proteins"
  if [ -d "$target_dir" ]; then
    mv "$file" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for file '$file'"
  fi
done
```

# FastANI for the 25 isolates

First we need to activate the gtdb environment:

```
conda activate gtdbtk-2.3.2
```

We create a new directory in the Results called FastANI:

```
mkdir /home/nik_arapitsas/Documents/Bacillus_project/Results/FastANI
```

We create a bash file containing the paths to the assemblies of each isolate. I ran this command while I was in the "/mnt/assemblies_repository" directory. We need a tab separated txt file that will have 2 columns: the first column will have the path to the FASTA file and the second column the six first letters of the filename (the isolate ID). We need to be in the "/mnt/assemblies_repository" to run this command:

```
find "$(pwd)" -type f -name "*.fasta" | awk -F'/' '{file=$NF; print $0 "\t" substr(file, 1, 6)}' > /home/nik_arapitsas/Documents/Bacillus_project/Results/FastANI/genome_list.txt
```

*In the "/mnt/assemblies_repository" directory there were more assemblies than the 25 (SRL307, SRL376 and SRL550). I manually deleted from the bash file the SRL307 and SRL376 isolates. I kept SRL550 for cross-validation. 

Go to the FastANI directory and run the FastANI using the GTDBTK database:

```
gtdbtk ani_rep --batchfile /home/nik_arapitsas/Documents/Bacillus_project/Results/FastANI/genome_list.txt --out_dir . --cpus 20
```

It ran in about 15 minutes.   


```
mkdir /home/nik_arapitsas/Documents/Bacillus_project/Results/FastANI_noresult_repeat
```

```
gtdbtk ani_rep --batchfile /home/nik_arapitsas/Documents/Bacillus_project/Results/FastANI/genome_list_no_result_isolates.txt --out_dir . --cpus 20
```

# AntiSMASH output assessment and visualisation

## Unzip all the .zip files containing the Antismash output

In the directory "/home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash" I had all the "assembly_Antismash".zip files and they were named x__assembly_Antismash.zip where x is the three digit number of the isolate. I wanted to extract the content of each zip file in separate folders in this directory and the folder shoud have the name SRLx_assembly_Antismash (SRL followed by the three digit number of each name x) and the zip files should have been removed.I created a script to automate the procedure: 

```
cd scripts/
touch unzip_antismash_files.sh
```

I put the code provided below to the script file:

```
for zipfile in *_assembly_Antismash.zip; do
    num=$(echo "$zipfile" | grep -oE '^[0-9]+')
    outdir="SRL${num}_assembly_Antismash"

    mkdir -p "$outdir"
    unzip "$zipfile" -d "$outdir"
    rm "$zipfile"
done
```

Made it executable by writing: 

```
chmod +x unzip_antismash_files.sh
```

I moved in the "/home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash" directory and run the script as shown below: 

```
cd /home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash

/home/nik_arapitsas/Documents/Bacillus_project/scripts/unzip_antismash_files.zip
```

## Count the antiSMASH regions for every isolate

I created the script count_antismash_regions.sh to go into every isolate folder and count the files that have the word "region" and are of .gbk type. Then, I went to the Antismash_Output directory and ran the script:  

```
cd Antismash_Output

/home/nik_arapitsas/Documents/Bacillus_project/scripts/count_antismash_regions.sh 
```

The next steps were performed in Rstudio for visualization. The file antismash_regions_graphs.R contains the code. 

## Movement of the antiSMASH data to meet the new server directory architecture 

```
cd /home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash/Antismash_Output

for dir in */; do
  dir="${dir%/}"
  prefix="${dir:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}"
  if [ -d "$target_dir" ]; then
    mv "$dir" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for folder '$dir'"
  fi
done
```

# Run quast in every assembly

I created a script named "quast_toeveryisolate.sh". 

I activated the quast environment: 

```
conda activate quast
```

I went to the folder where the assemblies are saved:

```
cd /media/sarlab/DATA/Bacillus_project/Assemblies/
```
I ran the script:

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/quast_toeveryisolate.sh
```

The path of the quast directories changed after moving the assemblies to meet the new server directory architechture that is described on the first section of this file.  

# Re-run unicycler in selected isolates with high contig number

## Isolate SRL662

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL662_assembly/SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_1.fq.gz -2 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL662_assembly/SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_2.fq.gz -l /media/sarlab/DATA/Bacillus_project/Assemblies/SRL662_assembly/SRL662_raw_data/A01_long.fastq -o /media/sarlab/DATA/Bacillus_project/Assemblies/SRL662_assembly/SRL662_new_assembly
--threads 8
```    

### Run quast for SRL662

```
conda activate quast
```
```
quast SRL662_new_assembly/assembly.fasta -o SRL662_new_assembly/SRL662_new_assembly_quast
``` 

## Isolate SRL368

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL368_assembly/SRL368_raw_data/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1.fq.gz -2 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL368_assembly/SRL368_raw_data/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2.fq.gz -l /media/sarlab/DATA/Bacillus_project/Assemblies/SRL368_assembly/SRL368_raw_data/368_bam.fastq -o /media/sarlab/DATA/Bacillus_project/Assemblies/SRL368_assembly/SRL368_new_assembly --threads 8
```    

### Run quast for SRL368

```
conda activate quast
```
```
quast SRL368_new_assembly/assembly.fasta -o SRL368_new_assembly/SRL368_new_assembly_quast
``` 

## Isolate SRL543

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL543_assembly/SRL543_raw_data/A08_FDSW210370234-1r_HLG2FDSX2_L1_1.fq.gz -2 /media/sarlab/DATA/Bacillus_project/Assemblies/SRL543_assembly/SRL543_raw_data/A08_FDSW210370234-1r_HLG2FDSX2_L1_2.fq.gz -l /media/sarlab/DATA/Bacillus_project/Assemblies/SRL543_assembly/SRL543_raw_data/A08_long.fastq -o /media/sarlab/DATA/Bacillus_project/Assemblies/SRL543_assembly/SRL543_new_assembly --threads 8
```    

### Run quast for SRL368

```
conda activate quast
```
```
quast SRL543_new_assembly/assembly.fasta -o SRL543_new_assembly/SRL543_new_assembly_quast
``` 

# Re-organization of the files (All the changes collected together)

In the directory "/media/sarlab/DATA/Bacillus_project" there were two directories called "proteins_Bacillus_project" and "proteins_Bacillus_th_israel". I moved them to the "DATA" directory temporaly to finish the formating of the main isolate directories in the "Bacillus_project" directory. 

Then I did the formating automatically: 

```
cd /media/sarlab/DATA/Bacillus_project
```

With the below code every folder will keep only the first six characters, which is the name of every isolate:

```
for d in */; do mv "$d" "${d:0:6}"; done
```

Then I used the code below to create a directory called "*Isolate code*_assembly" for every isolate and then move every assembly-associated file in it: 

```
for d in */; do 
  mkdir "$d"/"${d:0:6}"_assembly
  find "$d" -mindepth 1 -maxdepth 1 -type f -exec mv {} "$d"/"${d:0:6}"_assembly/ \;
done
```

Move all the proteins files in the respect directory:

```
cd /home/nik_arapitsas/Desktop/Test/proteins_Bacillus_project

for file in *; do
  prefix="${file:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}/${prefix}_proteins"
  if [ -d "$target_dir" ]; then
    mv "$file" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for file '$file'"
  fi
done
```

```
cd /home/nik_arapitsas/Downloads/Test

for file in *; do
  prefix="${file:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}/${prefix}_proteins"
  if [ -d "$target_dir" ]; then
    mv "$file" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for file '$file'"
  fi
done
```

## Transfer of antiSMASH files

```
cd /home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash/Antismash_Output

for dir in */; do
  dir="${dir%/}"
  prefix="${dir:0:6}"
  target_dir="/media/sarlab/DATA/Bacillus_project/${prefix}"
  if [ -d "$target_dir" ]; then
    mv "$dir" "$target_dir/"
  else
    echo "Warning: Target directory '$target_dir' not found for folder '$dir'"
  fi
done
```

# Trying a different assembly approach with flye and then Unicycler

## Isolate SRL662

Unicycler number of contigs: 35

### Inspect the Long reads

#### LongQC

```
conda activate perfect_assembly
conda install bioconda::longqc
```
```
cd media/sarlab/DATA/Bacillus_project/SRL662
```
```
python /home/labguest/software/LongQC/longQC.py sampleqc -x pb-sequel -o ./SRL662_longqc ./SRL662_raw_data/output.bc1011_1--bc1011_1.subreads.bam
```

#### NanoPlot

```
conda activate nanoplot
NanoPlot -o SRL662_nanoplot -p SRL662 --title "SRL662 NanoPlot Analysis" --ubam ./SRL662_raw_data/output.bc1011_1--bc1011_1.subreads.bam
```

#### Filtering (again)

Convert bam to fq: 

```
samtools fastq ./SRL662_raw_data/output.bc1011_1--bc1011_1.subreads.bam > ./SRL662_raw_data/A01_long_unfiltered.fastq
```
```
filtlong --min_length 1000 --keep_percent 95 ./SRL662_raw_data/A01_long_unfiltered.fastq > ./SRL662_raw_data/A01_long_filtered.fastq
```

##### Try filtering with the help of Illumina reads

```
filtlong -1 ./SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 ./SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --min_length 1000 --keep_percent 90 --target_bases 300000000 --trim --split 500 ./SRL662_raw_data/A01_long_unfiltered.fastq > ./SRL662_raw_data/A01_long_illumina_filtered.fastq
```

**Installation of the seqkit tool in the perfect_assembly environment**

```
conda install bioconda::seqkit
```

```
seqkit stats ./SRL662_raw_data/A01_long_unfiltered.fastq 
seqkit stats ./SRL662_raw_data/A01_long_filtered.fastq 
```

#### Run unicycler with unfiltered and filtered reads

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long_filtered.fastq -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_unicycler_filtered --threads 20
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long_unfiltered.fastq -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_unicycler_unfiltered --threads 20
flye --pacbio-raw ./SRL662_raw_data/A01_long_filtered.fastq --out-dir ./SRL662_flye_filtered_reads --genome-size 4.2m --threads 20
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l ./SRL662_raw_data/A01_long_filtered.fastq --existing_long_read_assembly ./SRL662_flye_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler --threads 20
```

Flye with the Illumina filtered long reads 

```
flye --pacbio-raw ./SRL662_raw_data/A01_long_illumina_filtered.fastq --out-dir ./SRL662_flye_illumina_filtered_reads --genome-size 4.2m --threads 20
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l ./SRL662_raw_data/A01_long_illumina_filtered.fastq --existing_long_read_assembly ./SRL662_flye_illumina_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_illumina_filtered_reads_unicycler --threads 20
```

### FastQC on the raw short reads

```
conda activate fastqc
cd /media/sarlab/DATA/Bacillus_project/SRL662/
mkdir SRL662_fastqc
fastqc ../SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_1.fq.gz ../SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastqc
```

### Fastp on the raw short reads

```
conda activate perfect_assembly
```
```
fastp -i /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_1.fq.gz -I /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_FDSW210370227-1r_HLG2FDSX2_L1_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -O /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --report_title "SRL662 fastp report" --unpaired1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_unpaired.fq.gz --unpaired2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_unpaired.fq.gz
```
### FastQC on the trimmed short reads

```
cd ..
mkdir SRL662_fastqc_trimmed
cd SRL662_fastqc_trimmed/
conda activate fastqc
fastqc ../SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz ../SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastqc_trimmed
```

### Run Flye on the long reads with the default parameters

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_20250520 --threads 15
```

I got 14 contigs.

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/
quast assembly.fasta -o ./SRL662_flye_assembly_20250520_quast
```

### Feed the Flye assembly with default parameters, to the Unicycler 

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520 --threads 16
```  

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/
quast assembly.fasta -o ./SRL662_flye_hybrid_assembly_20250520_quast
```

### Feed the Flye assembly with default parameters, to the Unicycler using --mode bold

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227
-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW21037
0227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_lo
ng.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20
250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_bold_2025
0602 --threads 23 --mode bold
```  

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_bold_20250602
quast assembly.fasta -o ./SRL662_flye_hybrid_assembly_bold_20250520_quast
```

### Try Flye with the --asm-coverage parameter and --meta parameter

**1) "--asm-coverage 50"**

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL662/
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_20250526 --genome-size 4.3m --asm-coverage 50 --threads 20
```

I got 204 contigs so I need to try a different command:

**2) "--meta"**

```
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_meta_20250526 --genome-size 4.3m --meta --threads 20
```

I got 22 contigs.

**3) "--asm-coverage 40"**

```
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_coverage40_20250526 --genome-size 4.3m --asm-coverage 40   --threads 20
```

I got 128 contigs.

**4) "--meta" without setting an estimated genome size**

```
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_meta_nogenomesize_20250526 --meta --threads 23
```
I got 22 contigs.

### Try unicycler on flye assembly with --asm-coverage 40

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_coverage40_20250526/
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_coverage40_assembly_20250530 --threads 23
```

### Use the Flye assembly that was produced with the "meta" parameter (case 2 from above) to the Unicycler

I will check if the 22 contigs will get reduced by using unicycler.

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_meta_20250526/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_meta_hybrid_assembly_20250526 --threads 23
```  

No! I got again 30 contigs as in the first case. So **for SRL662, flye with the default settings is the best option until now**.

### Try to run a second unicycler on SRL662 using the first unicycler assembly as input

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_unicycleroutput_hybrid_assembly_20250528 --threads 23
```

### Try polishing SRL662 with pilon 

```
bwa index /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/assembly.fasta
```

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/
```
``` 
bwa mem -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz | samtools sort -o SRL662_flye_hybrid_assembly_20250520_illumina_aligned.bam
```

```
samtools index SRL662_flye_hybrid_assembly_20250520_illumina_aligned.bam
```

```
pilon -Xmx200G --genome /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/assembly.fasta --frags SRL662_flye_hybrid_assembly_20250520_illumina_aligned.bam --changes --output SRL662_flye_hybrid_assembly_20250520_polished_first_round --outdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/
```

```
wc -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_pilon/SRL662_flye_hybrid_assembly_20250520_polished_first_round.changes
```

No changes have happened! I will try also racon.

### Try polishing SRL662 with racon

Racon uses the long reads for polishing, while pilon uses the short reads.

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_racon
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_racon
cp /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_racon
```

```
minimap2 -ax map-pb -t 23 assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq | samtools sort -o SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.sam
```
```
samtools view -S -b SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.sam > SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.bam
```
```
samtools index SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.bam
```
```
samtools view -h SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.bam > SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.sam
```

```
racon -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq SRL662_flye_hybrid_assembly_20250520_pacbio_aligned.sam assembly.fasta > SRL662_flye_hybrid_assembly_20250520_racon_polished_1st_round.fasta
```

No change on the contig number! 

### Try to use racon on the flye assembly before feeding it in the unicycler

#### 1st round

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_racon
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_racon/
minimap2 -ax map-pb -t 23 assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > SRL662_flye_assembly_20250520_pacbio_aligned.sam
racon -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq SRL662_flye_assembly_20250520_pacbio_aligned.sam assembly.fasta > SRL662_flye_assembly_20250520_racon1.fasta
```

#### 2nd round

```
minimap2 -ax map-pb -t 23 SRL662_flye_assembly_20250520_racon1.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > SRL662_flye_assembly_20250520_racon1_pacbio_aligned.sam
racon -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq SRL662_flye_assembly_20250520_racon1_pacbio_aligned.sam SRL662_flye_assembly_20250520_racon1.fasta > SRL662_flye_assembly_20250520_racon2.fasta
```

#### 3rd round

```
minimap2 -ax map-pb -t 23 SRL662_flye_assembly_20250520_racon2.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > SRL662_flye_assembly_20250520_racon2_pacbio_aligned.sam
racon -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq SRL662_flye_assembly_20250520_racon2_pacbio_aligned.sam SRL662_flye_assembly_20250520_racon2.fasta > SRL662_flye_assembly_20250520_racon3.fasta
```

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_racon/
conda activate quast
quast SRL662_flye_assembly_20250520_racon3.fasta -o ./SRL662_flye_assembly_20250520_racon3_quast
```

### Run unicycler using the racon polished flye assembly

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_3xracon_hybrid_assembly_20250528
```

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly SRL662_flye_assembly_20250520_racon3.fasta  -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_3xracon_hybrid_assembly_20250528 --threads 23
```

### Try polishing the flye assembly using pilon before using it in the unicycler

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/
cp ../assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/  
```

```
bwa index /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/assembly.fasta
```
``` 
bwa mem -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz | samtools sort -o SRL662_flye_assembly_20250520_illumina_aligned.bam
```

```
samtools index SRL662_flye_assembly_20250520_illumina_aligned.bam
```

```
pilon -Xmx200G --genome /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/assembly.fasta --frags SRL662_flye_assembly_20250520_illumina_aligned.bam --changes --output SRL662_flye_assembly_20250520_pilon1 --outdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/
```
```
wc -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_pilon/SRL662_flye_assembly_20250520_pilon1.changes
```
There were no changes. I stoped trying to reduce the contigs more. 

### Try to estimate genome size of SRL662 with jellyfish

Activate the environment where jellyfish is located:

```
conda activate genome_estimation
```

Create a directory where the jellyfish output will be located:

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_jellyfish
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_jellyfish
```

Move the Illumina short read raw data in this folder for making the code easier:

```
cp /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz .
cp /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz .
```

Combine the Illumina reads in one fq.gz file:

```
cat A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz > SRL662_illumina_reads_combined.fq.gz
```

Unzip the fq.gz file to obtain an .fq file:

```
gunzip -c SRL662_illumina_reads_combined.fq.gz > SRL662_illumina_reads_combined.fq
```

Run jellyfish for different k sizes:

**k = 15**

```
jellyfish count -C -m 15 -s 100000000000 -t 20 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k15.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k15.jf > SRL662_k15_reads.histo
```

**k = 17**

```
jellyfish count -C -m 17 -s 100M -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k17.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k17.jf > SRL662_k17_reads.histo
```

**k = 19**

```
jellyfish count -C -m 19 -s 100000000 -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k19.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k19.jf > SRL662_k19_reads.histo
```

**k = 21**

```
jellyfish count -C -m 21 -s 100M -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k21.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k21.jf > SRL662_k21_reads.histo
```

**k = 23**

```
jellyfish count -C -m 23 -s 100M -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k23.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k23.jf > SRL662_k23_reads.histo
```

**k = 15**

```
jellyfish count -C -m 15 -s 100M -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k15_lower_s.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k15_lower_s.jf > SRL662_k15_reads_lower_s.histo
```

**k = 13**

```
jellyfish count -C -m 13 -s 100M -t 23 SRL662_illumina_reads_combined.fq -o SRL662_mer_counts_k13.jf
```
```
jellyfish histo -t 20 SRL662_mer_counts_k13.jf > SRL662_k13_reads.histo
```

The k = 15 gave the best model fit and the estimated genome was 4,224,919 bp. 

### Run flye using in the parameters the estimated genome length for SRL662

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_estgenomesize_20250530 --genome-size 4224919 --threads 23
```

I got 14 contigs.

You can try --plasmids + You can filter reads to achieve certain coverage!

### Try reducing min overlap value

```
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_estgenomesize_minoverlap2000_20250530 --min-overlap 2000 --genome-size 4224919 --threads 23
```

I got 15 contigs. 

### Try to subsample long reads

#### For 100x coverage

```
mkdir ./SRL662_subsampled_100x
cd ./SRL662_subsampled_100x/  
```
```
seqtk sample -s100 ../SRL662_raw_data/A01_long.fastq 0.073 > SRL662_long_subsampled_100x.fastq
```
```
flye --pacbio-raw SRL662_long_subsampled_100x.fastq --out-dir ../SRL662_flye_assembly_estgenomesize_subsampled_100x_20250530 --genome-size 4224919 --threads 23
```

**I got 58 contigs**

##### Quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_100x_20250530
quast assembly.fasta -o ./SRL662_flye_assembly_estgenomesize_subsampled_100x_20250530_quast
```

#### For 860x coverage

```
mkdir ./SRL662_subsampled_860x
cd ./SRL662_subsampled_860x/  
```

```
seqtk sample -s100 ../SRL662_raw_data/A01_long.fastq 0.6278 > SRL662_long_subsampled_860x.fastq
```
```
flye --pacbio-raw SRL662_long_subsampled_860x.fastq --out-dir ../SRL662_flye_assembly_estgenomesize_subsampled_860x_20250601 --genome-size 4224919 --threads 23
```

**I got 12 contigs** 

#### For ~1000x coverage

```
mkdir ./SRL662_subsampled_1000x
cd ./SRL662_subsampled_1000x/  
```

```
seqtk sample -s100 ../SRL662_raw_data/A01_long.fastq 0.7447 > SRL662_long_subsampled_1000x.fastq
```
```
flye --pacbio-raw SRL662_long_subsampled_1000x.fastq --out-dir ../SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601 --genome-size 4224919 --threads 23
```

I got 12 contigs

##### Quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601
quast assembly.fasta -o ./SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601_quast
```

#### For ~600x coverage

```
mkdir ./SRL662_subsampled_600x
cd ./SRL662_subsampled_600x/  
```

```
seqtk sample -s100 ../SRL662_raw_data/A01_long.fastq 0.44682 > SRL662_long_subsampled_600x.fastq
```
```
flye --pacbio-raw SRL662_long_subsampled_600x.fastq --out-dir ../SRL662_flye_assembly_estgenomesize_subsampled_600x_20250601 --genome-size 4224919 --threads 23
```

I got 17 contigs

#### For ~400x coverage

```
mkdir ./SRL662_subsampled_400x
cd ./SRL662_subsampled_400x/  
```

```
seqtk sample -s100 ../SRL662_raw_data/A01_long.fastq 0.29788 > SRL662_long_subsampled_400x.fastq
```
```
flye --pacbio-raw SRL662_long_subsampled_400x.fastq --out-dir ../SRL662_flye_assembly_estgenomesize_subsampled_400x_20250601 --genome-size 4224919 --threads 23
```

I got 16 contigs

**The best is the x100 coverage asssembly (created using the x1000 coverage reads subsample). I will try this for unicycler**

#### Unicycler with the 100x coverage assembly

```
cd ..
mkdir ./SRL662_flye_hybrid_assembly_100xcoverage_20250601
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_20250601 --threads 23
```

#### Quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_20250601/
quast assembly.fasta -o ./SRL662_flye_hybrid_assembly_100xcoverage_20250601_quast
```

**It is wrong to use the whole long reads in the unicycler**

**I will run it again and use the subsampled reads.**

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_20250601 --threads 23
```

**Again 30 contigs**

### Try bbnorm for the short reads to achieve uniform coverage of 100x of every region

```
mkdir ../SRL662_norm_illumina_reads
cd ../SRL662_norm_illumina_reads
bbnorm.sh in1=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz in2=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz out1=./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm.fq.gz out2=./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm.fq.gz target=100 min=5 threads=23 
```
#### Run unicycler with the subsampled long-read assembly and using the normalized short reads

```
unicycler -1 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm.fq.gz -2 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads_20250602 --threads 23
```

### Try bbnorm for the short reads to achieve uniform coverage of 150x of every region

```
cd ../SRL662_norm_illumina_reads
bbnorm.sh in1=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz in2=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz out1=./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm150.fq.gz out2=./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm150.fq.gz target=150 mindepth=5 threads=23
```

I have not used them yet. I think maybe I need to use smaller coverage.

#### Run unicycler with the subsampled long-read assembly and using the normalized short reads x150

```
unicycler -1 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm150.fq.gz -2 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm150.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly_graph.gfa -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads150_20250602 --threads 23
```

I got 26 contigs from the 30! 


### Try bbnorm for the short reads to achieve uniform coverage of 50x of every region

```
cd ../SRL662_norm_illumina_reads
bbnorm.sh in1=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz in2=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz out1=./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm50.fq.gz out2=./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm50.fq.gz target=50 mindepth=5 threads=23
```

#### Run unicycler with the subsampled long-read assembly and using the normalized short reads x50

```
unicycler -1 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm50.fq.gz -2 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm50.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads50_20250602 --threads 23
```

I got 26 contigs from the 30! 

#### Quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads50_20250602
quast assembly.fasta -o ./SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads50_20250602_quast
```

#### Run unicycler with the subsampled long-read assembly and using the normalized short reads x50

```
unicycler -1 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm50.fq.gz -2 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm50.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly_graph.gfa -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads50gfa_20250602 --threads 23
```


### Try bbnorm for the short reads to achieve uniform coverage of 30x of every region

```
cd ../SRL662_norm_illumina_reads
bbnorm.sh in1=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz in2=/media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz out1=./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm30.fq.gz out2=./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm30.fq.gz target=30 mindepth=5 threads=23
```

#### Run unicycler with the subsampled long-read assembly and using the normalized short reads x30

```
unicycler -1 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed_norm30.fq.gz -2 ./A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed_norm30.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_subsampled_1000x/SRL662_long_subsampled_1000x.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_estgenomesize_subsampled_1000x_20250601/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_100xcoverage_normIlluminareads30_20250602 --threads 23
```

I got 43 contigs. 

### Try canu as an alternative long-read assembler to flye

```
canu -p SRL662 -d /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_canu_assembly genomeSize=4224919 -pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq
```

I got 464 contigs

### Try raven as an alternative assembler

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raven_assembly_20250603 
```

```
raven -t 18 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raven_assembly_20250603/SRL662_raven_assembly.fasta
```

**I got 9 contigs!**

#### Quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raven_assembly_20250603/
quast SRL662_raven_assembly.fasta -o ./SRL662_raven_assembly_20250603_quast
```

#### Try unicycler on the raven assembly

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raven_assembly_20250603/SRL662_raven_assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raven_hybrid_assembly_20250603 --threads 23
```

Again 30 contigs :(

### Try hybrid SPAdes

```
spades.py -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --pacbio /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_hybrid_spades_assembly_20250604 -t 23
```

Try the isolate mode:

```
spades.py -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --pacbio /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --isolate -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_hybrid_spades_isolate_assembly_20250604 -t 23
```

trusted assembly

```
spades.py -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --pacbio /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --isolate --trusted-contigs /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_20250520_racon/SRL662_flye_assembly_20250520_racon3.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_hybrid_spades_isolate_trc_assembly_20250604 -t 23
```

### Install inspector

```
conda create -n inspector \
  python=3.8 \
  pysam \
  statsmodels=0.10.1 \
  minimap2=2.15 \
  samtools=1.9 \
  flye=2.8.3 \
  -c bioconda -c conda-forge
conda activate inspector
```
```
conda install pandas=1.0.5
conda install numpy=1.21
```

```
git clone https://github.com/ChongLab/Inspector.git
```

#### Try inspector for SRL662

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_inspector_20250604
```

```
inspector.py -c /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/assembly.fasta -r /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_inspector_20250604 -t 23 --datatype clr
```

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_inspector_20250604
```

```
inspector.py -c /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/assembly.fasta -r /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_inspector_20250604 -t 23 --datatype clr
```

```
inspector-correct.py -i /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/SRL662_flye_assembly_inspector_20250604/ --datatype pacbio-raw -t 23
```

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly contig_corrected.fa -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_inspcorrected_20250604 --threads 23 
```

### Try minimap2/miniasm

```
 mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_minimap2_miniasm
```
```
minimap2 -x ava-pb -t 23 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq | gzip -1 > /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_minimap2_miniasm/SRL662_minimap2_overlaps.paf.gz
```

```
miniasm -f /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_minimap2_miniasm/SRL662_minimap2_overlaps.paf.gz > /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_minimap2_miniasm/SRL662_assembly.gfa
```

```
awk '/^S/{print ">" $2 "\n" $3}' SRL662_assembly.gfa | fold > SRL662_miniasm_contigs.fasta
```

#### Try the less strict miniasm

```
miniasm -c 2 -e 2 -f /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq SRL662_minimap2_overlaps.paf.gz > SRL662_assembly_less_strict.gfa
```
```
awk '/^S/{print ">" $2 "\n" $3}' SRL662_assembly_less_strict.gfa | fold > SRL662_miniasm_less_strict_contigs.fasta
```

```
conda activate quast
quast SRL662_miniasm_contigs.fasta -o ./quast
quast SRL662_miniasm_less_strict_contigs.fasta -o ./quast_less_strict
```

### Try nextdenovo

```
conda activate perfect_assembly
conda install bioconda::nextdenovo
mkdir SRL662_nextdenovo
cd SRL662_nextdenovo
```

NextDenovo needs the path of the reads' file to a .fofn file. We create that with the following command:

```
ls /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > SRL662_input.fofn
```

Then we need to create a .cfg file to feed it in the nexdenovo assembler. I got the code from this webpage: https://nextdenovo.readthedocs.io/en/latest/TEST1.html, and created the SRL662_nextdenovo_run.cfg script file. I changed the following parameters: 

* job_type = local  (from sge)

* read_type = clr (from ont)

* workdir = ./ (to have the output in my current directory)

* genome_size = 4m  

Then, I run the nextdenovo assembler with the following command:

```
nextDenovo SRL662_nextdenovo_run.cfg
```

### Try smartdenovo

```
conda install bioconda::smartdenovo
```
```
mkdir SRL662_smartdenovo
cd SRL662_smartdenovo 
smartdenovo.pl -p SRL662 -t 20 -c 1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq > SRL662.mak 
make -f SRL662.mak
```

I got 21 contigs

### Try WTDBG2

```
conda activate perfect_assembly
conda install bioconda::wtdbg
```

```
mkdir SRL662_wtdbg2
cd SRL662_wtdbg2
wtdbg2 -x sq -g 4.2m -t 20 -i /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq -fo SRL662  
wtpoa-cns -t 20 -i SRL662.ctg.lay.gz -fo SRL662.ctg.fa
```

#### Quast 

```
conda activate quast
quast SRL662.ctg.fa -o ./SRL662_wtdbg2_quast
cd SRL662_wtdbg2_quast/
vim report.txt
```

**225 contigs but the only long-read assembler that gave a good length: 4591404 bp**

#### Feed the WTDBG2 output in unicycler

```
mkdir SRL662_wtdbg2_unicycler
cd SRL662_wtdbg2_unicycler
```
```
conda activate perfect_assembly
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --existing_long_read_assembly ../SRL662_wtdbg2/SRL662.ctg.fa -o ./ --threads 23 
```

### Try autocycler

#### Install autocycler

```
conda activate perfect_assembly
conda install bioconda::autocycler 
```

#### Run autocycler

**1) Subsample** 

```
autocycler subsample --reads /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long.fastq --out_dir ./SRL662_autocycler/SRL662_subsampled_reads --genome_size 4224919
```

**2) Run assemblers for each subset**

i) Flye

```
mkdir SRL662_assemblies
```
```
for i in 01 02 03 04; do
        /home/nik_arapitsas/Documents/Bacillus_project/scripts/flye.sh SRL662_subsampled_reads/sample_"$i".fastq SRL662_assemblies/flye_"$i" 23 4224919 
done
```


## Isolate SRL368

Unicycler number of contigs: 15 (Christos), 18 (Nikos)

### Fastp on the raw short reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL368/
mkdir SRL368_fastp
cd SRL368_fastp
```
```
fastp -i /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1.fq.gz -I /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -O /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_trimmed.fq.gz --report_title "SRL368 fastp report" --unpaired1 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_unpaired.fq.gz --unpaired2 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_unpaired.fq.gz
```

### Run Flye on the long reads

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL368_raw_data/368_bam.fastq --out-dir ./SRL368_flye_assembly_20250520 --threads 15
```

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_assembly_20250520/
quast assembly.fasta -o ./SRL368_flye_assembly_20250520_quast
```

### Use the Flye assembly to the Unicycler

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/368_bam.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_assembly_20250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_hybrid_assembly_20250520 --threads 16
```  

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_hybrid_assembly_20250520/
quast assembly.fasta -o ./SRL368_flye_hybrid_assembly_20250520_quast
```

### Try Flye with the --asm-coverage 50 parameter

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL368/
flye --pacbio-raw ./SRL368_raw_data/368_bam.fastq --out-dir ./SRL368_flye_assembly_20250526 --genome-size 7m --asm-coverage 50 --threads 20
```

#### Filtering of the long reads

Convert bam to fq: 

```
samtools fastq sample_368_m64270e_230219_013217.subreads.bam > SRL368_long_unfiltered.fastq
```
```
filtlong --min_length 1000 --keep_percent 95 SRL368_long_unfiltered.fastq > SRL368_long_filtered.fastq
```

### Run flye and unicycler with the filtered long-reads

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq --out-dir /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_filtered_reads --threads 23
```

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq --out-dir /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_filtered_reads --threads 23
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_filtered_reads_unicycler --threads 23
```

### Try running flye with --plasmids parameter

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq --out-dir /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_filtered_reads_plasmids --plasmids --threads 18
```

### Try hybrid SPAdes

```
spades.py -1 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_trimmed.fq.gz --pacbio /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_hybrid_spades_assembly -t 18
```

```
spades.py -1 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_fastp/sample_368_FKDN230011102-1A_HYVFYDSX3_L3_2_trimmed.fq.gz --pacbio /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_raw_data/SRL368_long_filtered.fastq --plasmid -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_hybrid_plasmid_spades_assembly -t 18
```


## Isolate SRL543

Unicycler number of contigs: 6 (Christos), 9 (Nikos)

### Fastp on the raw short reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL543/
mkdir SRL543_fastp
cd SRL543_fastp
```
```
fastp -i /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_FDSW210370234-1r_HLG2FDSX2_L1_1.fq.gz -I /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_FDSW210370234-1r_HLG2FDSX2_L1_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -O /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_2_trimmed.fq.gz --report_title "SRL543 fastp report" --unpaired1 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_1_unpaired.fq.gz --unpaired2 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_2_unpaired.fq.gz
```

### Run Flye on the long reads

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL543_raw_data/A08_long.fastq --out-dir ./SRL543_flye_assembly_20250521 --threads 20
```
```
flye --pacbio-raw ./SRL543_raw_data/A08_long.fastq --out-dir ./SRL543_flye_assembly_20250521 --genome-size 5m --threads 20
```

There were no assembly on the previous cases and a message saying "ERROR: No disjointigs were assembled - please check if the read type and genome size parameters are correct
[2025-05-21 18:07:33] ERROR: Pipeline aborted".

I searched and found that the parameter **--asm-coverage 50** could help.

```
flye --pacbio-raw ./SRL543_raw_data/A08_long.fastq --out-dir ./SRL543_flye_assembly_20250525 --genome-size 5m --asm-coverage 50 --threads 20
```
**It worked!**

```
grep -c "^>" /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250525/assembly.fasta 
```
**It gave 4 contigs!**

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250525/
quast assembly.fasta -o ./SRL543_flye_assembly_20250525_quast
```

### Use the Flye assembly to the Unicycler

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250525/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_hybrid_assembly_20250525 --threads 20
```  

```
grep -c "^>" /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_hybrid_assembly_20250525/assembly.fasta 
```

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_hybrid_assembly_20250525/
quast assembly.fasta -o ./SRL543_flye_hybrid_assembly_20250525_quast
```

**It has 1 contig and now is ok**

### Filter the PacBio raw reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data
```

Convert bam to fq: 

```
samtools fastq output.1069_1--1069_1.subreads.bam > A08_long_unfiltered.fastq
```

Run filtlong:

```
filtlong --min_length 1000 --keep_percent 95 A08_long_unfiltered.fastq > ./A08_long_filtered.fastq
```

### Run flye and unicycler with the filtered long-reads

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_long_filtered.fastq --genome-size 5m --asm-coverage 50 --out-dir /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_filtered_reads --threads 23
```

**I got 1 contig only by flye this time.**

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_long_filtered.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_filtered_reads_unicycler --threads 23
```

## Isolate SRL389

Unicycler number of contigs: 15 (Christos), 14 (Nikos)

### Fastp on the raw short reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL389/
mkdir SRL389_fastp
cd SRL389_fastp
```
```
fastp -i /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1.fq.gz -I /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -O /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2_trimmed.fq.gz --report_title "SRL389 fastp report" --unpaired1 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1_unpaired.fq.gz --unpaired2 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2_unpaired.fq.gz
```

### Re-run unicycler

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/389_bam.fastq -o /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_assembly_20250527 --threads 23
```

I got 14 contigs.

### Run Flye on the long reads

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL389_raw_data/389_bam.fastq --out-dir ./SRL389_flye_assembly_20250527 --threads 23
```

I got 3 contigs! 

### Run unicycler using the Flye assembly

Now I will do the unicycler using the Flye assembly:

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/389_bam.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_assembly_20250527/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_hybrid_assembly_20250527 --threads 23
```

I got 8 contigs.

### Run quast

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_hybrid_assembly_20250527/
quast assembly.fasta -o ./SRL389_flye_hybrid_assembly_20250527_quast
```

```
rfplasmid --species Bacillus --input /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_hybrid_assembly_20250527/ --threads 23 --out /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_hybrid_assembly_20250527/SRL389_flye_hybrid_assembly_20250527_rfplasmid
```

### Filter the PacBio raw reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data
```

Convert bam to fq: 

```
samtools fastq sample_389_m64270e_230219_013217.subreads.bam > SRL389_long_unfiltered.fastq
```

Run filtlong:

```
filtlong --min_length 1000 --keep_percent 95 SRL389_long_unfiltered.fastq > ./SRL389_long_filtered.fastq
```

### Run flye and unicycler with the filtered long-reads

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/SRL389_long_filtered.fastq --out-dir /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads --threads 23
```

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_fastp/sample_389_FKDN230011106-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_raw_data/SRL389_long_filtered.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads_unicycler --threads 23
```

## Isolate SRL340

Unicycler number of contigs: 5 (Christos), 6 (Nikos)

### Fastp on the raw short reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL340/
mkdir SRL340_raw_data
mkdir SRL340_fastp
cd SRL340_fastp
```
```
fastp -i /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1.fq.gz -I /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2.fq.gz -o /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -O /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2_trimmed.fq.gz --report_title "SRL340 fastp report" --unpaired1 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1_unpaired.fq.gz --unpaired2 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2_unpaired.fq.gz
```

### Re-run unicycler

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/340_bam.fastq -o /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_assembly_20250527 --threads 23
```

I got 6 contigs.

### Run Flye on the long reads

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL340_raw_data/340_bam.fastq --out-dir ./SRL340_flye_assembly_20250527 --threads 23
```

I got 11 contigs.

### Try Flye with different minimum overlap

```
flye --pacbio-raw ./SRL340_raw_data/340_bam.fastq --out-dir ./SRL340_flye_assembly_minoverlap_10k_20250527 --threads 23 --min-overlap 10000
```

It gave 90 contigs. So I will proceed with the default. 

### Run unicycler using the Flye assembly

Now I will do the unicycler using the Flye assembly:

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/340_bam.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_assembly_20250527/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_hybrid_assembly_20250527 --threads 23
```

I got 6 contigs.

### Run quast

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_hybrid_assembly_20250527/
quast assembly.fasta -o ./SRL340_flye_hybrid_assembly_20250527_quast
```

### Filter the PacBio raw reads

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data
```

Convert bam to fq: 

```
samtools fastq sample_340_m64270e_230219_013217.subreads.bam > SRL340_long_unfiltered.fastq
```

Run filtlong:

```
filtlong --min_length 1000 --keep_percent 95 SRL340_long_unfiltered.fastq > ./SRL340_long_filtered.fastq
```

### Run flye and unicycler with the filtered long-reads

```
flye --pacbio-raw /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/SRL340_long_filtered.fastq --out-dir /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_filtered_reads --threads 23
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_fastp/sample_340_FKDN230011111-1A_HYVFYDSX3_L3_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_raw_data/SRL340_long_filtered.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_filtered_reads/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_flye_filtered_reads_unicycler --threads 23
```

# Search for plasmids

## Install rfplasmid

```
conda create --name plasmid_search
conda activate plasmid_search
conda install bioconda::rfplasmid
rfplasmid --initialize
```

## Install plasmidfinder

```
conda activate plasmid_search
conda install bioconda::plasmidfinder
download-db.sh
```

Upon its installation, plasmidfinder had the not "Please run download-db.sh to download the PlasmidFinder database to /opt/miniconda3/envs/plasmid_search/share/plasmidfinder-2.1.6/database.
If you have a database in custom path, please use plasmidfinder.py with the option -p." in the end. So that is why I ran the command download-db.sh. 

**SOS** The database to run plasmidfinder is located to this path: 

```
/opt/miniconda3/envs/plasmid_search/share/plasmidfinder-2.1.6/database/
```

Remember to create the output folder as well before running the plasmdifinder or it will print "Input Error: Output dirctory does not exist!".

## Use rfplasmid and plasmidfinder for the isolates

### SRL368 

```
rfplasmid --species Bacillus --input /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly_20230327/ --threads 23 --out /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly_20230327/SRL368_assembly_20230327_rfplasmid
```


```
mkdir /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly_20230327/SRL368_assembly_20230327_plasmidfinder 
plasmidfinder.py -i /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly_20230327/SRL368_assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly_20230327/SRL368_assembly_20230327_plasmidfinder -p /opt/miniconda3/envs/plasmid_search/share/plasmidfinder-2.1.6/database/ 
```

### SRL662

```
rfplasmid --species Bacillus --input /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/ --threads 23 --out /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_rfplasmid
```

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_plasmidfinder 
plasmidfinder.py -i /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_hybrid_assembly_20250520/SRL662_flye_hybrid_assembly_20250520_plasmidfinder -p /opt/miniconda3/envs/plasmid_search/share/plasmidfinder-2.1.6/database/ 
```

### SRL224

```
rfplasmid --species Bacillus --input /media/sarlab/DATA/Bacillus_project/SRL224/SRL224_assembly --threads 23 --out /media/sarlab/DATA/Bacillus_project/SRL224/SRL224_assembly_rfplasmid
```

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL224/SRL224_assembly_plasmidfinder 
plasmidfinder.py -i /media/sarlab/DATA/Bacillus_project/SRL224/SRL224_assembly/SRL224_assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL224/SRL224_assembly_plasmidfinder -p /opt/miniconda3/envs/plasmid_search/share/plasmidfinder-2.1.6/database/ 
```


## Use ragtag

Install ragtag:

```
conda install -c bioconda ragtag
```

Run ragtag correct: 

```
cd /media/sarlab/DATA/Bacillus_project/SRL662 
mkdir SRL662_ragtag
cd SRL662_ragtag
mkdir SRL662_ragtag_correct
cd SRL662_ragtag_correct
```
```
ragtag.py correct /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta 
```

Run ragtag scaffold:

```
cd ../..
mkdir SRL662_ragtag_scaffold
cd SRL662_ragtag_scaffold
```
```
ragtag.py scaffold /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta 
```

Run ragtag patch:

```
cd ../..
mkdir SRL662_ragtag_patch
cd SRL662_ragtag_patch
```

```
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta
```

I did not run this command. **I preffered to use the ragtag correct output for ragtag scaffold and then run ragtag patch**: 

```
cd ..
mkdir SRL662_ragtag_correct_scaffold
cd SRL662_ragtag_correct_scaffold
```

Use the corrected assembly for the scaffold command: 

```
ragtag.py scaffold /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct/ragtag_output/ragtag.correct.fasta 
```

Then run ragtag patch using the ragtag correct output:

```
cd ..
mkdir SRL662_ragtag_patch_usingcorrected
cd SRL662_ragtag_patch_usingcorrected
```

```
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct/ragtag_output/ragtag.correct.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta  
```

Try inverting the position of SRL658 and SRL662 assemblies in the command:

```
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct/ragtag_output/ragtag.correct.fasta
```

Quast to compare SRL658 with SRL662 assembly:

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_invert/ragtag_output/SRL662vsSRL658_quast
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_invert/ragtag_output/SRL662vsSRL658_quast
```
```
quast /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_invert/ragtag_output/ragtag.patch.fasta -r /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta -o ./
```

Try unicycler using the ragtag.patch.fasta file as "existing long-read assembly":

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag
mkdir SRL662_ragtag_patch_usingcorrected_invert_unicycler
cd SRL662_ragtag_patch_usingcorrected_invert_unicycler
```

```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fastp/A01_FDSW210370227-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long_filtered.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_invert/ragtag_output/ragtag.patch.fasta -o /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_invert_unicycler --threads 20
```

Using scaffold as input for target in ragtap patch:

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag
mkdir SRL662_ragtag_patch_usingcorrected_scaffold
cd SRL662_ragtag_patch_usingcorrected_scaffold
```

```
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta
```
```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_scaffold/ragtag_output/SRL662_ragtag_patch_usingcorrected_scaffold_patch2
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_scaffold/ragtag_output/SRL662_ragtag_patch_usingcorrected_scaffold_patch2
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_scaffold/ragtag_output/ragtag.patch.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta
```

Try inverted

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag
mkdir SRL662_ragtag_patch_usingcorrected_scaffold_inverted
cd SRL662_ragtag_patch_usingcorrected_scaffold_inverted
```

```
ragtag.py patch /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta
```

```
nucmer -p SRL662_asm_comparison /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_patch_usingcorrected_scaffold_inverted/ragtag_output/ragtag.patch.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta
dnadiff -p SRL662_dnadiff_out -d SRL662_asm_comparison.delta
```

```
nucmer -p SRL658vsSRL656_asm_comparison /media/sarlab/DATA/Bacillus_project/SRL656/SRL656_assembly/SRL656_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta
dnadiff -p SRL658vsSRL656_dnadiff_out -d SRL658vsSRL656_asm_comparison.delta
```
```
nucmer -p SRL662last_asm_comparison /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta
dnadiff -p SRL662last_dnadiff_out -d SRL662last_asm_comparison.delta
```

**The SRL662, SRL658 and SRL656 are probably the same isolate. In the case of SRL662 the final assembly is more fragmented and the total length appears smaller, maybe due to repeats that interfered with the complete assembly of SRL662** 


## Run TGS-Gapcloser

```
mkdir /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_tgsgapcloser
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_tgsgapcloser
```

Convert the long-reads from .fastq to .fasta:

```
seqtk seq -a /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_raw_data/A01_long_filtered.fastq > A01_long_filtered.fasta
```

Run tgsgapcloser with racon error correction:

```
tgsgapcloser --scaff /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta --tgstype pb --reads A01_long_filtered.fasta --output /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_tgsgapcloser --racon /opt/miniconda3/envs/perfect_assembly/bin/racon --thread 20 
```

Run tgsgapcloser without error correction:

```
mkdir SRL662_tgsgapcloser_noerrorcorr
cd SRL662_tgsgapcloser_noerrorcorr
```

```
tgsgapcloser --scaff /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta --tgstype pb --reads ../A01_long_filtered.fasta --output . --ne --thread 20 >pipe.log 2>pipe.err
```

## Try using GapCloser

```
conda activate perfect_assembly
conda install bioconda::soapdenovo2-gapcloser
```

```
GapCloser -a /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta -b /media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta -o SRL662_ragtag_correct_scaffold_gapcloser -t 20
```

## Try using Fgap

```
conda create -n fgap -c conda-forge -c bioconda fgap
```

```
FGAP "-d /home/nik_arapitsas/Desktop/SRL662_temporary/SRL662_ragtag/SRL662_ragtag_correct_scaffold/ragtag_output/ragtag.scaffold.fasta -a '/media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta' -o ./ -t 18"
```
```
FGAP "-d /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fgap/_2.fasta -a '/media/sarlab/DATA/Bacillus_project/SRL658/SRL658_assembly/SRL658_assembly.fasta' -o ./ -t 18"
```
```
FGAP "-d /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_fgap/_2.fasta -a '/media/sarlab/DATA/Bacillus_project/SRL656/SRL656_assembly/SRL656_assembly.fasta' -o ./ -t 18"
```

# Antismash

## Install antiSMASH

```
conda create -n antismash -c conda-forge -c bioconda antismash
conda activate antismash
download-antismash-databases
```

## Run antiSMASH

```
antismash --genefinding-tool prodigal --taxon bacteria --cpus 18 --output-dir ./ /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta
```

```
antismash --genefinding-tool prodigal --taxon bacteria --cpus 18 --output-dir ./SRL662_antismash_test /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta --fullhmmer --clusterhmmer --tigrfam --asf --cb-general --cc-mibig --cb-subclusters --cb-knownclusters --pfam2go --rre --smcog-trees            
```

**I used the webpage to upload each assembly and download the online output: https://antismash.secondarymetabolites.org/#!/start** 


# BUSCO on the assemblies

I created the script *mkdir_for_busco.sh* to create a directory for the busco output in each isolate folder, and to run busco for every assembly. 


# Run prodigal for isolates SRL389, SRL543 and SRL662

## SRL389

```
conda activate perfect_assembly
cd /media/sarlab/DATA/Bacillus_project/SRL389
mkdir SRL389_proteins
cd SRL389_proteins
prodigal -i /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads_unicycler/assembly.fasta -o SRL389_gene_coordinates.gff -a SRL389_proteins.faa -d SRL389_genes.fna
```

## SRL543

```
cd /media/sarlab/DATA/Bacillus_project/SRL543
mkdir SRL543_proteins
cd SRL543_proteins
prodigal -i /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_filtered_reads_unicycler/assembly.fasta -o SRL543_gene_coordinates.gff -a SRL543_proteins.faa -d SRL543_genes.fna
```

## SRL662

```
cd /media/sarlab/DATA/Bacillus_project/SRL662
mkdir SRL662_proteins
cd SRL662_proteins
prodigal -i /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta -o SRL662_gene_coordinates.gff -a SRL662_proteins.faa -d SRL662_genes.fna
```

# Orthofinder with the new assemblieas and adjustments

## Re-run Orthofinder with every isolate

```
conda activate orthofinder
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_proteins -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder
```

It took about 15 minutes to complete. It produced 741.5 Mb of output files. 

## Run Orthofinder while keeping out isolates that are the same species and are isolated from the exact same sample

**Groups of same species from the same sample:**

Sample MIES02:	SRL215,  SRL218,  SRL224  

Sample MIES03:	SRL221,  SRL244

Sample S10:	SRL398,   SRL342

Sample MTR:	SRL656,    SRL658,    SRL662

* MIES02:

grep -c "^>" SRL215_proteins.faa 
6135

grep -c "^>" SRL218_proteins.faa 
6097

grep -c "^>" SRL224_proteins.faa 
6097

**I selected SRL215**

* MIES03:

grep -c "^>" SRL221_proteins.faa 
3813

grep -c "^>" SRL244_proteins.faa 
3813

**I selected SRL244**


* S10:

grep -c "^>" SRL342_proteins.faa 
6619

grep -c "^>" SRL398_proteins.faa 
6591

**I selected BOTH (pretty interesting)**

SRL342 (S10c1 -> NA1/2) and SRL398 (S10a5b1 -> NA) both belong to the same species. The only difference is that they were cultivated in different media for their consecutive isolations and re-inoculations. That could indicate media-associated diversification that can have an impact in the genome. 

* MTR:

grep -c "^>" SRL656_proteins.faa 
4260

grep -c "^>" SRL658_proteins.faa 
4260

grep -c "^>" SRL662_proteins.faa 
3811

As these all are probably the same, and as SRL662 has lower number of proteins and its assembly is fragmented, I decided to select one from the other two. **I just picked SRL656.**

So, I got rid of SRL221, SRL662, SRL658, SRL218 and SRL224 protein files from the "Bacillus_project_proteins" directory and moved them in a directory called "Bacillus_project_redundant_proteins". 

```
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_proteins -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder_filtered
```

## Create Graphs

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs
```

1) **Number of Isolate-Specific Orthogroups per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/species_specific_orthogroups.txt  
```

The graph was designed in Rstudio. The script is located in the following path:

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/orthofinder_graphs.R
```

2) **Percentage of genes from each isolate assigned to orthogroups**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

The same sorting with the list above will be used to provide plots that could be easily comparable:

First, extract the species order from the previous output file with species-specific orthogroups:

```
awk '{print $1}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/species_specific_orthogroups.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/isolates.txt 
```

Then use it to assign the percentage of genes in orthofroups in the desirable order: 

```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/percofgenes_inogs_per_isolate_unsorted.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/isolates.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/percofgenes_inogs_per_isolate.txt
```

3) **Genes with orthogroups in all or any isolates**

With the code below when counting the partially shared orthogroups we do not count the core orthogroups.

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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/orthogroupcount_in_isolates.txt
```

With the code below when counting the partially shared orthogroups we count the core orthogroups as well. This is better for creating a bar plot where the bars of all and any orthogroups will be overlapping. 

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    species[i] = substr($i, 1, 6); # Keep only first 6 letters of species name
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      if (core) core_count[i]++;  # Count genes in core orthogroups
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break} # Count genes in shared orthogroups, including core
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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/orthogroupcount_in_isolates.txt
```

# Orthofinder between isolates of the same species that are originated from the same samples

## Run Orthofinder for SRL215 SRL218 SRL224

```
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/B_thurigiensis -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/B_thurigiensis_orthofinder
```

## Run Orthofinder for SRL342 and SRL398

```
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder
```

nucmer -p SRL342_vs_SRL398_comparison /media/sarlab/DATA/Bacillus_project/SRL342/SRL342_assembly/SRL342_assembly.fasta /media/sarlab/DATA/Bacillus_project/SRL398/SRL398_assembly/SRL398_assembly.fasta

dnadiff -p SRL342_vs_SRL398_comparison_dnadiff_out -d SRL342_vs_SRL398_comparison.delta


# Use roary to make pangenome comparisons between isolates

```
conda create -n roary
conda activate roary
conda install -c conda-forge -c bioconda roary
```

```
roary -e -n -v -f ./SRL342_vs_SRL398_roary *.gff -p 20
```

### 1) Get the species-specific Orthogroups for SRL398 from the Orthofinder of 25 isolates

```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignoring last column
{
  target = 20;  # Column for isolate T (adjust if needed)
  if ($target > 0) {  # Ensure isolate T has genes
    is_species_specific = 1;

    # Check if any other isolate (excluding column 6) has genes
    for (i=2; i<=NF-2; i++) {  # NF-2 to ignore the last column
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups.txt
```

### 2) Get the species-specific Genes for SRL398 by searching using species-specific Orthogroups

First I need to use tab as delimiter in the txt file: 

```
sed 's/ \+/	/g' SRL398_species_specific_orthogroups.txt > SRL398_species_specific_orthogroups_with_tabs.txt
```

In order to get the isolate-specific genes I used the following command:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from species_specific_orthogroups.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $20;  # Extract the gene names from the fifth column
        print orthogroup, genes;  # Print orthogroup name and corresponding gene names
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/SRL398_specific_ogs_and_genes/SRL398_species_specific_genes.txt
```

Get as output the genes in a list:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from SRL_398_species_specific_orthogroups_with_tabs.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $20;  # Extract the gene names from the fifth column
        split(genes, gene_array, ",");  # Split the gene names into an array
        for (i in gene_array) {
            print orthogroup, gene_array[i];  # Print orthogroup and gene in two columns
        }
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/SRL398_specific_ogs_and_genes/SRL398_species_specific_genes.txt
```

Extract the protein sequence of each isolate specific gene: 

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/SRL398_extract_isolate_specific_genes.sh
```

### 1) Get the species-specific Orthogroups for SRL398 from the Orthofinder of just SRL342 and SRL398 

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes
```

```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignoring last column
{
  target = 3;  # Column for SRL398 (third column)
  if ($target > 0) {  # Ensure SRL398 has genes
    is_species_specific = 1;

    # Check all other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 to skip "Total" column
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups.txt
```
```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignore last column
{
  target = 2;  # Column for SRL342
  if ($target > 0) {  # Check target has genes
    is_species_specific = 1;

    # Check other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 because last column is "Total"
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_orthogroups.txt
```

### 2) Get the species-specific Genes for SRL398 and SRL342 by searching using species-specific Orthogroups

First I need to use tab as delimiter in the txt file: 

```
sed 's/ \+/	/g' SRL398_species_specific_orthogroups.txt > SRL398_species_specific_orthogroups_with_tabs.txt
sed 's/ \+/	/g' SRL342_species_specific_orthogroups.txt > SRL342_species_specific_orthogroups_with_tabs.txt
```

In order to get the isolate-specific genes I used the following command:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from species_specific_orthogroups.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $3;  # Extract the gene names from the fifth column
        print orthogroup, genes;  # Print orthogroup name and corresponding gene names
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_genes.txt
```
```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from species_specific_orthogroups.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $2;  # Extract the gene names from the fifth column
        print orthogroup, genes;  # Print orthogroup name and corresponding gene names
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_genes.txt
```

Get as output the genes in a list:

```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from SRL_398_species_specific_orthogroups_with_tabs.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $3;  # Extract the gene names from the fifth column
        split(genes, gene_array, ",");  # Split the gene names into an array
        for (i in gene_array) {
            print orthogroup, gene_array[i];  # Print orthogroup and gene in two columns
        }
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_genes.txt
```
```
awk -F'\t' '
NR==FNR {species_specific[$1]; next}  # Read orthogroups from SRL_340_species_specific_orthogroups_with_tabs.txt into an array
{
    orthogroup = $1;  # Get orthogroup name from the first column of Orthogroups.tsv
    if (orthogroup in species_specific) {  # If orthogroup exists in species-specific list
        genes = $2;  # Extract the gene names from the fifth column
        split(genes, gene_array, ",");  # Split the gene names into an array
        for (i in gene_array) {
            print orthogroup, gene_array[i];  # Print orthogroup and gene in two columns
        }
    }
}
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_orthogroups_with_tabs.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/Orthogroups/Orthogroups.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_redundant_proteins/Paenibacillus_sp_orthofinder/Results_Jun26/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_genes.txt
```

Extract the protein sequence of each isolate specific gene: 

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/SRL398_pair_extract_isolate_specific_genes.sh
/home/nik_arapitsas/Documents/Bacillus_project/scripts/SRL342_pair_extract_isolate_specific_genes.sh
```

# Run Orthofinder for the prodigal output

```
conda activate orthofinder
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_proteins -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder
```

## Create Graphs

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs
```

1) **Number of Isolate-Specific Orthogroups per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/species_specific_orthogroups.txt  
```

The graph was designed in Rstudio. The script is located in the following path:

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/orthofinder_graphs.R
```

2) **Percentage of genes from each isolate assigned to orthogroups**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

The same sorting with the list above will be used to provide plots that could be easily comparable:

First, extract the species order from the previous output file with species-specific orthogroups:

```
awk '{print $1}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/species_specific_orthogroups.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/isolates.txt 
```

Then use it to assign the percentage of genes in orthofroups in the desirable order: 

```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/percofgenes_inogs_per_isolate_unsorted.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/isolates.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/percofgenes_inogs_per_isolate.txt
```

3) **Genes with orthogroups in all or any isolates**

With the code below when counting the partially shared orthogroups we do not count the core orthogroups.

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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/orthogroupcount_in_isolates.txt
```

With the code below when counting the partially shared orthogroups we count the core orthogroups as well. This is better for creating a bar plot where the bars of all and any orthogroups will be overlapping. 

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    species[i] = substr($i, 1, 6); # Keep only first 6 letters of species name
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      if (core) core_count[i]++;  # Count genes in core orthogroups
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break} # Count genes in shared orthogroups, including core
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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Graphs/orthogroupcount_in_isolates.txt
```

## Find the unique orthogroups for SRL342 and SRL398

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/SRL342_SRL398_specific_ogs_and_genes
```

```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignoring last column
{
  target = 14;  # Column for SRL342 (fourteenth column)
  if ($target > 0) {  # Ensure SRL398 has genes
    is_species_specific = 1;

    # Check all other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 to skip "Total" column
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups.txt
```
```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignore last column
{
  target = 20;  # Column for SRL342
  if ($target > 0) {  # Check target has genes
    is_species_specific = 1;

    # Check other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 because last column is "Total"
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul01/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups.txt
```

# Annotate the genomes using Bakta

## Install Bakta

```
conda create -n bakta
conda install -c conda-forge -c bioconda bakta
bakta_db download --output  --type full
```

** The databse has an unziped size of 89.6 Gb

## Run Bakta

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/run_bakta.sh
```
```
bakta --db /media/sarlab/DATA/Bacillus_project/SRL662/db --prefix SRL662_bakta --output /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_bakta --threads 20 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/assembly.fasta
```

# Copy for the assemblies SRL340, SRL368, SRL389 and SRL662 the contigs larger than 1000 bp in separate files

## SRL340 

```
cd /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_assembly
seqkit seq -m 1000 -g SRL340_assembly.fasta > SRL340_assembly_filtered.fasta 
seqkit stats SRL340_assembly_filtered.fasta 
```

## SRL368 

```
cd /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly
seqkit seq -m 1000 -g SRL368_assembly.fasta > SRL368_assembly_filtered.fasta 
seqkit stats SRL368_assembly_filtered.fasta 
```

## SRL389

```
cd /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads_unicycler
seqkit seq -m 1000 -g assembly.fasta > SRL389_assembly_filtered.fasta 
seqkit stats SRL389_assembly_filtered.fasta 
```

## SRL662 

```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler
seqkit seq -m 1000 -g assembly.fasta > SRL662_assembly_filtered.fasta 
seqkit stats SRL662_assembly_filtered.fasta 
```

# Run Bakta on the selected contigs for these four isolates

## SRL662

```
bakta --db /media/sarlab/DATA/databases/bakta_v6.0/db --prefix SRL662_bakta --output /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_filtered_bakta --threads 20 /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_filtered_reads_unicycler/SRL662_assembly_filtered.fasta
```

## The rest of isolates 

```
bakta --db /media/sarlab/DATA/databases/bakta_v6.0/db --prefix SRL340_bakta --output /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_filtered_bakta --threads 20 /media/sarlab/DATA/Bacillus_project/SRL340/SRL340_assembly/SRL340_assembly_filtered.fasta

bakta --db /media/sarlab/DATA/databases/bakta_v6.0/db --prefix SRL368_bakta --output /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_filtered_bakta --threads 20 /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_assembly/SRL368_assembly_filtered.fasta

bakta --db /media/sarlab/DATA/databases/bakta_v6.0/db --prefix SRL389_bakta --output /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_filtered_bakta --threads 20 /media/sarlab/DATA/Bacillus_project/SRL389/SRL389_flye_filtered_reads_unicycler/SRL389_assembly_filtered.fasta
```

# Run Orthofinder for the bakta output 

## Re-run Orthofinder with every isolate

```
conda activate orthofinder
orthofinder -f /media/sarlab/DATA/Bacillus_project/Bacillus_project_proteins -t 20 -o /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder
```

## Create Graphs

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs
```

1) **Number of Isolate-Specific Orthogroups per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==9 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/species_specific_orthogroups.txt  
```

The graph was designed in Rstudio. The script is located in the following path:

```
/home/nik_arapitsas/Documents/Bacillus_project/scripts/orthofinder_graphs.R
```

2) **Percentage of genes from each isolate assigned to orthogroups**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==5 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/percofgenes_inogs_per_isolate_unsorted.txt  
```

The same sorting with the list above will be used to provide plots that could be easily comparable:

First, extract the species order from the previous output file with species-specific orthogroups:

```
awk '{print $1}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/species_specific_orthogroups.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/isolates.txt 
```

Then use it to assign the percentage of genes in orthofroups in the desirable order: 

```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/percofgenes_inogs_per_isolate_unsorted.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/isolates.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/percofgenes_inogs_per_isolate.txt
```

3) **Genes with orthogroups in all or any isolates**

With the code below when counting the partially shared orthogroups we do not count the core orthogroups.

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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/orthogroupcount_in_isolates.txt
```

With the code below when counting the partially shared orthogroups we count the core orthogroups as well. This is better for creating a bar plot where the bars of all and any orthogroups will be overlapping. 

```
awk '
NR==1 {
  for (i=2; i<=NF-1; i++) {
    species[i] = substr($i, 1, 6); # Keep only first 6 letters of species name
    core_count[i] = 0;
    shared[i] = 0;
  }
  next
}
{
  core=1;
  for (i=2; i<=NF-1; i++) if ($i==0) core=0; # Check if this is a core orthogroup

  for (i=2; i<=NF-1; i++) {
    if ($i > 0) {
      if (core) core_count[i]++;  # Count genes in core orthogroups
      for (j=2; j<=NF-1; j++) {
        if (j != i && $j > 0) {shared[i]++; break} # Count genes in shared orthogroups, including core
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
' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/orthogroupcount_in_isolates.txt
```

4) **Number of Unassigned Genes per Isolate - Genetic Novelty**

```
awk -F'\t' 'NR==1 {for(i=2; i<=NF; i++) species[i]=substr($i, 1, 6)} NR==6 {for(i=2; i<=NF; i++) print species[i], $i}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sort -k2,2n > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/unasigned_genes_perisolate.txt  
```
```
awk 'NR==FNR{a[$1]=$2; next} $1 in a {print $1, a[$1]}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/unasigned_genes_perisolate.txt /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/isolates.txt > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Graphs/unasigned_genes_perisolate_sorted.txt
```

## Find the unique orthogroups of SRL342 and SRL398 

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/SRL342_SRL398_specific_ogs_and_genes
```

```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignoring last column
{
  target = 14;  # Column for SRL342 (fourteenth column)
  if ($target > 0) {  # Ensure SRL398 has genes
    is_species_specific = 1;

    # Check all other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 to skip "Total" column
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/SRL342_SRL398_specific_ogs_and_genes/SRL342_species_specific_orthogroups.txt
```
```
awk -F'\t' '
NR==1 {for(i=2; i<=NF-1; i++) species[i]=$i; next}  # Store species names, ignore last column
{
  target = 20;  # Column for SRL398
  if ($target > 0) {  # Check target has genes
    is_species_specific = 1;

    # Check other species columns except target
    for (i=2; i<=NF-1; i++) {  # NF-1 because last column is "Total"
      if (i != target && $i > 0) {
        is_species_specific = 0;
        break;
      }
    }

    if (is_species_specific) {
      print $1, species[target];  # Print orthogroup ID & species name
    }
  }
}' /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/Orthogroups/Orthogroups.GeneCount.tsv > /media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jul03/SRL342_SRL398_specific_ogs_and_genes/SRL398_species_specific_orthogroups.txt
```

# Use PanACoTA for Pangenome analysis

```
mkdir /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota
cd /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota
mkdir SRL368_panacota
mkdir SRL179_panacota
mkdir SRL337_panacota
mkdir SRL543_panacota
```

## Install ncbi-genome-download package to get the closest relative sequences from the NCBI

```
conda create -n ncbi_genome_download
conda activate ncbi_genome_download
conda install -c bioconda ncbi-genome-download
```

## Run ncbi-genome-download command 

i) SRL368

```
ncbi-genome-download bacteria --assembly-accessions GCF_020775355.1,GCF_022663675.1,GCF_000021305.1,GCF_002559825.1,GCF_001182785.1 --formats fasta --flat-output --output-folder /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL368_panacota
gunzip *.fna.gz
```

ii) SRL179

```
ncbi-genome-download bacteria --assembly-accessions GCF_030123405.1,GCF_000612625.1,GCF_023715505.1 --formats fasta --flat-output --output-folder /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL179_panacota
gunzip /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL179_panacota/*.fna.gz
```

iii) SRL337

```
ncbi-genome-download bacteria --assembly-accessions GCF_003581585.1,GCF_003581585.1,GCF_008180415.1 --formats fasta --flat-output --output-folder /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL337_panacota
gunzip /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL337_panacota/*.fna.gz
```

iv) SRL543

```
ncbi-genome-download bacteria --assembly-accessions GCF_024706405.1,GCF_020171945.1,GCF_000473245.1 --formats fasta --flat-output --output-folder /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL543_panacota
gunzip /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL543_panacota/*.fna.gz
```

There was a newer version for the assembly GCF_024706405.1, the GCF_024706405.2

```
ncbi-genome-download bacteria --assembly-accessions GCF_024706405.2 --formats fasta --flat-output --output-folder /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL543_panacota
gunzip /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL543_panacota/*.fna.gz
```

## PANACOTA FOR SRL368

**I) ASSEMBLE**

```
cd /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL368_panacota
```
* PaNACoTA needs a list of file names as input. It can be created with the following script:

```
for f in *.fna *.fasta; do g="${f%.*}"; ext="${f##*.}"; clean=$(echo "$g" | tr '.' '_'); newname="${clean}.${ext}"; mv "$f" "$newname"; echo "$newname"; done > SRL368_list_genomes.lst
```

```
PanACoTA annotate -d /media/sarlab/DATA/Bacillus_project/Bacillus_project_panacota/SRL368_panacota -r ./annotation_output -l SRL368_list_genomes.lst -n BATH --threads 20
```






