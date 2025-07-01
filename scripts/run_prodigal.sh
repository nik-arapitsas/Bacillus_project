#!/bin/bash

conda activate perfect_assembly

base_dir="/media/sarlab/DATA/Bacillus_project"

for dir in "${base_dir}"/SRL662/; do
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
    
    # Create prodigal directory
    prodigal_dir="${dir}${shortname}_prodigal"
    mkdir -p "$prodigal_dir"
    
    echo "Found assembly: $assembly_file"
    echo "Created BUSCO dir: $prodigal_dir"

    # Run busco
    
    busco -i "$assembly_file" -o "$(basename "$prodigal_dir")" --out_path "$(dirname "$prodigal_dir")" -l firmicutes_odb10 -m genome -c 18 -f
    
    prodigal -i "$assembly_file" -o "${dir}${shortname}_gene_coordinates.gff" -a "${dir}${shortname}_proteins.faa" -d "${dir}${shortname}_genes.fna"

done

