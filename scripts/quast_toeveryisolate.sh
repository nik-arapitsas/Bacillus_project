#!/bin/bash -l

# Name: quast_toeveryisolate.sh
# Purpuse: perform quast to multiple isolates
# Author: Nikolaos Arapitsas
# Date: 07/05/2025
# Usage: ./quast_toeveryisolate.sh assemblies_directories.txt 

conda activate quast

for dir in *_assembly/; do
    if [ -d "$dir" ]; then
        fasta_file="$dir/*_assembly.fasta"
        output_dir="$dir/${dir:0:6}_quast"

        if [ -f "$fasta_file" ]; then
            echo "Running quast on $fasta_file..."
            quast "$fasta_file" -o "$output_dir"
        else
            echo "Warning: $fasta_file not found, skipping."
        fi
    fi
done

