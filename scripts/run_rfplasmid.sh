#!/bin/bash
######################################################################################################
# script name: run_rfplasmid.sh
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to run the rfplasmid command for every assembly in the 
# /media/sarlab/DATA/Bacillus_project directory
######################################################################################################
# !!! MAKE SURE TO USE: conda activate rfplasmid BEFORE RUNNING THE SCRIPT !!!
# usage:./run_rfplasmid.sh
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/run_rfplasmid.sh
######################################################################################################

base_dir="/media/sarlab/DATA/Bacillus_project"

for dir in "${base_dir}"/SRL*/; do
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
    
    rfplasmid_dir="${dir}${shortname}_assembly_rfplasmid"
    
    echo "Found assembly: $assembly_file"

    # Run rfplasmid
    
    rfplasmid --species Bacillus --input "$assembly_dir" --threads 23 --out "$rfplasmid_dir"

done

