#!/bin/bash

# Usage: ./panacota_extract_unique_accessory_genes.sh SRL179

X=$1

# Create output directory
mkdir -p ${X}_accessory

# Locate required files
lstfile=$(find ./${X}_pangenome -type f -name "PanGenome*.tsv.lst")
mapfile=$(find ./${X}_annotation_output -type f -name "LSTINFO-*list_genomes.lst")

# Loop through each genome ID used internally by PanACoTA
cut -f1 "$mapfile" | while read -r internal_id; do

  # Get the original filename (second column)
  orig_name=$(awk -v id="$internal_id" '$1 == id {print $2}' "$mapfile")

  # Determine output-friendly isolate ID
  if [[ $orig_name == SRL* ]]; then
    outname=$(echo "$orig_name" | cut -d'_' -f1)
  elif [[ $orig_name == GCF* ]]; then
    outname=$(echo "$orig_name" | awk -F'_' '{print $1 "_" $2}')
  else
    outname=$orig_name
  fi

  # Extract unique genes present only in this genome
  awk -v target="$internal_id" -v OFS="\t" '
  {
    keep = 1
    genes = ""
    for (i = 2; i <= NF; i++) {
      if ($i != "-") {
        if (index($i, target) != 1) {
          keep = 0
          break
        } else {
          genes = genes ? genes OFS $i : $i
        }
      }
    }
    if (keep && genes != "")
      print $1, genes
  }' "$lstfile" > ${X}_accessory/${outname}_unique_genes.txt

rm -f ${X}_accessory/orig_name_unique_genes.txt

done
