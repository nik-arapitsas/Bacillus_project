#!/bin/bash
######################################################################################################
# script name: collect_bakta_proteins.sh
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to collect all the 25 proteome files from bakta output in the  
# /media/sarlab/DATA/Bacillus_project_proteins directory
######################################################################################################
# usage:./collect_bakta_proteins.sh
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/collect_bakta_proteins.sh
######################################################################################################


base_dir="/media/sarlab/DATA/Bacillus_project"
protein_dir="${base_dir}/Bacillus_project_proteins"
mkdir -p "$protein_dir"

# Loop through the isolates (assuming they are in SRL662 or other folders directly inside base_dir)
for bakta_dir in "${base_dir}"/SRL*/*_bakta; do
    if [[ -d "$bakta_dir" ]]; then
        shortname=$(basename "$bakta_dir" _bakta)
        protein_file="${bakta_dir}/${shortname}_bakta.faa"
        if [[ -f "$protein_file" ]]; then
            cp "$protein_file" "${protein_dir}/"
            echo "Copied $protein_file to $protein_dir"
        else
            echo "Warning: $protein_file not found, skipping."
        fi
    fi
done