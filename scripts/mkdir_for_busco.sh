
for dir in /home/nik_arapitsas/Desktop/test/SRL*/; do
    
    shortname=$(basename "$dir")

    # Enter the directory
    cd "$dir" || { echo "Failed to enter $dir"; continue; }
    
    # Create the busco folder
    busco_dir="${shortname}_busco"
    mkdir -p "$busco_dir" && echo "Created $busco_dir in $dir"
    
    # Return to original directory (optional, but good practice)
    cd - >/dev/null

done


   