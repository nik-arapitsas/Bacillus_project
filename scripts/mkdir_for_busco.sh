
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
    
    assembly_file=$(ls "${assembly_dir}"/*assembly.fasta 2>/dev/null | head -1)
    
    [[ -f "$assembly_file" ]] || { echo "Error: No assembly file in $assembly_dir"; continue; }
    
    # Create BUSCO directory
    busco_dir="${dir}${shortname}_busco"
    mkdir -p "$busco_dir"
    
    echo "Found assembly: $assembly_file"
    echo "Created BUSCO dir: $busco_dir"

    # Run busco
    
    busco -i "$assembly_file" -o "$(basename "$busco_dir")" --out_path "$(dirname "$busco_dir")" -l firmicutes_odb10 -m genome -c 18 \
    
done

