#!/bin/bash

# Set directories
source_base="/media/sarlab/DATA/Bacillus_project"
target_dir="/media/sarlab/DATA/Bacillus_project/Bacillus_project_proteins"

# Create target directory
mkdir -p "$target_dir"

# Process each SRL directory
for srl_dir in "${source_base}"/SRL[0-9][0-9][0-9]*; do
    # Extract base SRL number (first 3 digits after SRL)
    srl_number=$(basename "$srl_dir" | grep -o 'SRL[0-9]\{3\}' | grep -o '[0-9]\{3\}')
    
    # Find the proteins file (handles various patterns)
    source_file=$(find "$srl_dir" -maxdepth 2 -name "SRL${srl_number}*proteins.faa" -type f -print -quit)
    
    if [[ -f "$source_file" ]]; then
        # Create standardized target filename
        target_file="${target_dir}/SRL${srl_number}_proteins.faa"
        cp -v "$source_file" "$target_file" && \
        echo "Copied: $(basename "$source_file") → SRL${srl_number}_proteins.faa"
    else
        echo "Warning: No proteins.faa file found in $srl_dir"
        # Debugging line (uncomment if needed):
        # find "$srl_dir" -name "*proteins.faa" -ls
    fi
done

echo "===================================="
echo "All protein files copied to: $target_dir"
echo "Original files remain in their locations"