#!/bin/bash

# Input files
fasta_file="/mnt/assemblies_repository/proteins_Bacillus_th_israel/SRL368_IMG_8028368599_proteins.faa"
genes_list="/home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/species_specific_genes_list.txt"
output_file="/home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder_SRL368_relatives/Results_Feb12/SRL368_specific_ogs_and_genes/filtered_protein_sequences.txt"

# Create an empty output file
> "$output_file"

# Normalize gene names: remove extra spaces and hidden characters
awk '{gsub(/\r/,""); gsub(/[ \t]+$/, ""); print $2}' "$genes_list" > gene_names.tmp

# Extract gene sequences
awk -v gene_list="gene_names.tmp" '
BEGIN {
    while ((getline < gene_list) > 0) {
        genes[$1] = 1;  # Store gene names in array
    }
    close(gene_list);
}
{
    if ($0 ~ /^>/) {  
        split($0, header, "#");  # Split at #, keep first part
        split(header[1], clean_name, " ");  # Remove trailing spaces
        gene_name = substr(clean_name[1], 2);  # Remove ">" character
        capturing = (gene_name in genes) ? 1 : 0;  # Check if gene is in list
        if (capturing) {
            if (sequence != "") print gene "\t" sequence;
            gene = gene_name;
            sequence = "";
        }
    } else if (capturing) {
        sequence = sequence $0;
    }
}
END {
    if (sequence != "") print gene "\t" sequence;
}
' "$fasta_file" >> "$output_file"

# Clean up
rm gene_names.tmp

echo "Extraction completed. Output saved to: $output_file"