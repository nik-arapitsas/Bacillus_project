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