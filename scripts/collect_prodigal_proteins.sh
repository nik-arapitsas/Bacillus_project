#!/bin/bash

base_dir="/media/sarlab/DATA/Bacillus_project"
protein_dir="${base_dir}/Bacillus_project_proteins"
mkdir -p "$protein_dir"

# Loop through the isolates (assuming they are in SRL662 or other folders directly inside base_dir)
for prodigal_dir in "${base_dir}"/SRL*/*_prodigal; do
    if [[ -d "$prodigal_dir" ]]; then
        shortname=$(basename "$prodigal_dir" _prodigal)
        protein_file="${prodigal_dir}/${shortname}_proteins.faa"
        if [[ -f "$protein_file" ]]; then
            cp "$protein_file" "${protein_dir}/"
            echo "Copied $protein_file to $protein_dir"
        else
            echo "Warning: $protein_file not found, skipping."
        fi
    fi
done