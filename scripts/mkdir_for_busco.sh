
#!/bin/bash

conda activate busco

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
    
    assembly_file="${assembly_dir}/assembly.fasta"
    
    if [[ ! -f "$assembly_file" ]]; then
        echo "Warning: assembly.fasta not found in $assembly_dir"
        continue
    fi
    
    # Create BUSCO directory
    busco_dir="${dir}${shortname}_busco"
    mkdir -p "$busco_dir"
    
    echo "Found assembly: $assembly_file"
    echo "Created BUSCO dir: $busco_dir"
    
    # Run busco
    
    busco -i "$assembly_file" -o "$busco_dir" -l bacteria_odb10 -m genome -c 8

done