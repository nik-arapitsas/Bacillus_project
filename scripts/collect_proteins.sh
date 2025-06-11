#!/bin/bash

# Set directories
source_base="/home/nik_arapitsas/Desktop/test/SRL368"
target_dir="/home/nik_arapitsas/Desktop/test/Bacillus_project_proteins"

# Create target directory
mkdir -p "$target_dir"

# Process each SRL directory
for srl_dir in "${source_base}"/SRL[0-9][0-9][0-9]; do
    # Extract SRL number (last 3 digits)
    srl_number=$(basename "$srl_dir" | grep -o '[0-9]\{3\}$')
    
    # Define file paths
    source_file="${srl_dir}/SRL${srl_number}_proteins/SRL${srl_number}_proteins.faa"
    target_file="${target_dir}/SRL${srl_number}_proteins.faa"
    
    # Verify and copy
    if [[ -f "$source_file" ]]; then
        cp -v "$source_file" "$target_file" && \
        echo "Copied: SRL${srl_number}_proteins.faa"
    else
        echo "Warning: Not found - ${srl_dir}/SRL${srl_number}_proteins/SRL${srl_number}_proteins.faa"
    fi
done

echo "===================================="
echo "All protein files copied to:"
echo "$target_dir"
echo "Original files remain in their locations"