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

### Run Flye on the long reads

```
cd ..
conda activate perfect_assembly
flye --pacbio-raw ./SRL662_raw_data/A01_long.fastq --out-dir ./SRL662_flye_assembly_20250520 --threads 15
```

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL662/SRL662_flye_assembly_20250520/
quast assembly.fasta -o ./SRL662_flye_assembly_20250520_quast
```

### Use the Flye assembly to the Unicycler

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

```
grep -c "^>" /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250521/assembly.fasta 
```

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250521/
quast assembly.fasta -o ./SRL543_flye_assembly_20250521_quast
```

### Use the Flye assembly to the Unicycler

```
conda activate perfect_assembly
```
```
unicycler -1 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_1_trimmed.fq.gz -2 /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_fastp/A08_FDSW210370234-1r_HLG2FDSX2_L1_2_trimmed.fq.gz -l /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_raw_data/A08_long.fastq --existing_long_read_assembly /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_assembly_20250520/assembly.fasta -o /media/sarlab/DATA/Bacillus_project/SRL368/SRL368_flye_hybrid_assembly_20250520 --threads 16
```  

#### Run quast on the assembly

```
conda activate quast 
```
```
cd /media/sarlab/DATA/Bacillus_project/SRL543/SRL543_flye_assembly_20250521
quast assembly.fasta -o ./SRL543_flye_hybrid_assembly_20250521_quast
```

