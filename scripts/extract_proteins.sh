temp_base="/home/nik_arapitsas/Downloads/Test/temp_folders"  # Define temp folder base path
output_folder="/mnt/assemblies_repository/proteins_Bacillus_project"  # Define where to save extracted files
mkdir -p "$temp_base" "$output_folder"  # Ensure required directories exist

for zipfile in /home/nik_arapitsas/Downloads/Test/SRL*.zip; do
    filename=$(basename "$zipfile")  # Extract only the filename
    prefix="${filename%%.*}"  # Extract prefix (SRL)

    temp_folder="$temp_base/$prefix"  # Set proper temp folder path
    mkdir -p "$temp_folder"  # Create temp folder

    unzip -d "$temp_folder" "$zipfile"  # Extract ZIP contents

    # Find the correct proteins.faa file, excluding cds_proteins.faa
    proteins_file=$(find "$temp_folder" -type f -name "*proteins.faa" ! -name "*cds_proteins.faa")

    if [[ -n "$proteins_file" ]]; then
        mv "$proteins_file" "$output_folder/${prefix}.faa"  # Rename and move extracted file
    else
        echo "Warning: No proteins.faa found in $zipfile"
    fi

    rm -rf "$temp_base"  # Delete temp folder
done
