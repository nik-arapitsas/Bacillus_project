#!/bin/bash

# Output file
output_file="region_counts.txt"

# Header (optional for R)
echo -e "isolate\tregion_count" > "$output_file"

# Loop through folders ending in assembly_Antismash
for dir in *assembly_Antismash/; do
    if [ -d "$dir" ]; then
        # Get count of .gbk files that contain "region"
        count=$(find "$dir" -type f -name "*region*.gbk" | wc -l)
        
        # Get first 6 letters of the folder name (remove trailing slash first)
        shortname=$(basename "$dir" | cut -c1-6)

        # Write to output file (tab-separated)
        echo -e "${shortname}\t${count}" >> "$output_file"
    fi
done

echo "Output written to $output_file"
