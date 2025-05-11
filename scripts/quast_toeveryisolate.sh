#!/bin/bash -l

# Name: quast_toeveryisolate.sh
# Purpuse: perform quast to multiple isolates
# Author: Nikolaos Arapitsas
# Date: 07/05/2025
# Usage: ./quast_toeveryisolate.sh assemblies_directories.txt 

for dir in *_assembly/; do
    if [ -d "$dir" ]; then
        basename="${dir%/}"
        output_dir="${dir}${basename:0:6}_quast"

        # Look for a file that ends with 'assembly.fasta' in the directory
        fasta_file=$(find "$dir" -maxdepth 1 -type f -name '*assembly.fasta' | head -n 1)

        if [ -n "$fasta_file" ]; then
            echo "Running quast on $fasta_file..."
            quast "$fasta_file" -o "$output_dir"
        else
            echo "Warning: No file ending in 'assembly.fasta' found in $dir, skipping."
        fi
    fi
done

