#!/bin/bash
######################################################################################################
# script name: run_bakta.sh
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to run the bakta annotation pipeline for every isolate in the 
# /media/sarlab/DATA/Bacillus_project directory
######################################################################################################
# !!! MAKE SURE TO USE: conda activate bakta BEFORE RUNNING THE SCRIPT !!!
# usage:./run_bakta.sh
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/run_bakta.sh
######################################################################################################

base_dir="/media/sarlab/DATA/Bacillus_project"

for dir in "${base_dir}"/SRL368/; do
    shortname=$(basename "$dir")
    echo "Processing $shortname"
    
    # Find the assembly directory (ending with either assembly or unicycler)
    assembly_dir=$(find "$dir" -maxdepth 1 -type d \( -name "*assembly" -o -name "*unicycler" \) -print -quit)
    
    if [[ -z "$assembly_dir" ]]; then
        echo "Warning: No assembly folder found in $shortname"
        continue
    fi
    
    assembly_file=$(ls "${assembly_dir}"/*assembly.fasta 2>/dev/null | head -1)
    
    [[ -f "$assembly_file" ]] || { echo "Error: No assembly file in $assembly_dir"; continue; }
    
    bakta_dir="${dir}${shortname}_bakta"
    
    echo "Found assembly: $assembly_file"
    echo "The output files will be found in the direcotry: $bakta_dir"
    echo "Their name will have the preffix: ${shortname}_bakta"

    # Run Bakta
    
    bakta --db /media/sarlab/DATA/Bacillus_project/SRL662/db --prefix "${shortname}_bakta" --output "${bakta_dir}" --threads 20 "$assembly_file"

done

