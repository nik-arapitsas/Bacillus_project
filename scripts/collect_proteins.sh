#!/bin/bash

# Set directories
base_dir="/home/nik_arapitsas/Desktop/test/SRL368"
target_dir="${base_dir}/Bacillus_project_proteins"

# Create target directory
mkdir -p "$target_dir"

# Process each SRL directory
for srl_dir in "${base_dir}"/SRL[0-9][0-9][0-9]/; do
    # Extract SRL number
    srl_number=$(basename "$srl_dir" | sed 's/SRL//')
    
    # Define file paths
    source_file="${srl_dir}SRL${srl_number}_proteins/SRL${srl_number}_proteins.faa"
    target_file="${target_dir}/SRL${srl_number}_proteins.faa"
    
    # Verify and copy
    if [[ -f "$source_file" ]]; then
        # Copy to central directory (preserves original in place)
        cp -v "$source_file" "$target_file"
    else
        echo "Warning: $source_file not found"
    fi
done

echo "All protein files copied to: $target_dir"
echo "Original files remain in their locations"