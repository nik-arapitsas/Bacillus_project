#!/bin/bash

# Output directory for cleaned files
mkdir -p cleaned_fastas

for f in GCF_*.fna SRL*; do
  [[ -e "$f" ]] || continue

  # Generate Anvi’o-friendly filename
  if [[ $f == GCF_* ]]; then
    newname=$(echo "$f" | sed -E 's/^(GCF_[^_]+)[^ ]*/\1/' | sed 's/[^a-zA-Z0-9_]/_/g').fna
  else
    newname=$(echo "$f" | sed -E 's/^(SRL[^_]+)_.*/\1.fasta/')
  fi
  echo "Renaming $f → $newname"
  mv "$f" "$newname"

  # Clean headers for Anvi’o (replace dots with underscores in first word)
  awk '/^>/ { 
          split($1, id, " "); 
          gsub(/\./, "_", id[1]); 
          print id[1]; 
        } 
       !/^>/ { print }' "$newname" > "cleaned_fastas/$newname"
done

echo "✅ Done. Cleaned FASTA files are in the 'cleaned_fastas/' directory."




# Define base output directories
profile_base="/media/sarlab/DATA/Bacillus_project/Bacillus_project_anvio/SRL179_anvio/SRL179_genomes_db_profiles"
summary_base="/media/sarlab/DATA/Bacillus_project/Bacillus_project_anvio/SRL179_anvio/SRL179_genomes_db_summaries"

# Make sure output directories exist
mkdir -p "$profile_base"
mkdir -p "$summary_base"

# Loop through each .db file
for db in *.db; do
  [ -e "$db" ] || continue  # Skip if no .db files
  prefix="${db%.*}"

  echo "🔄 Processing $prefix"

  # Step 1: Create blank profile
  anvi-profile -c "$db" \
               --blank-profile \
               --sample-name "$prefix" \
               -o "$profile_base/${prefix}_genome_profile" \
               -T 20

  # Step 2: Summarize the profile
  anvi-summarize -c "$db" \
                 -p "$profile_base/${prefix}_genome_profile/PROFILE.db" \
                 -o "$summary_base/${prefix}_genome_summary"

  echo "✅ Finished $prefix"
done

echo "🎉 All done! Profiles and summaries are in:"
echo " - Profiles: $profile_base"
echo " - Summaries: $summary_base"