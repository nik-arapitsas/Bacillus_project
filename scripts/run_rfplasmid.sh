#!/bin/bash

base_dir="/home/nik_arapitsas/Desktop/test"

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
    
    # Create rfplasmid directory
    rfplasmid_dir="${dir}${shortname}_assembly_rfplasmid"
    mkdir -p "$rfplasmid_dir"
    
    echo "Found assembly: $assembly_file"
    echo "Created rfplasmid dir: $rfplasmid_dir"

    # Run rfplasmid
    
    rfplasmid --species Bacillus --input "$assembly_dir" --threads 23 --out "$rfplasmid_dir"
    
    # Find the timestamped output subdirectory
    ts_dir=$(find "$rfplasmid_dir" -mindepth 1 -maxdepth 1 -type d -name "${shortname}_assembly_rfplasmid_*" | head -1)

    if [[ -n "$ts_dir" ]]; then
        echo "Moving contents from $ts_dir to $rfplasmid_dir"
        mv "$ts_dir"/* "$rfplasmid_dir"
        rmdir "$ts_dir"
    else
        echo "Warning: Timestamped directory not found in $rfplasmid_dir"
    fi

done

