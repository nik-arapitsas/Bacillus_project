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