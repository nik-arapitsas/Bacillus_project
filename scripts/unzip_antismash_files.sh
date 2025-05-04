for zipfile in *_assembly_Antismash.zip; do
    num=$(echo "$zipfile" | grep -oE '^[0-9]+')
    outdir="SRL${num}_assembly_Antismash"

    mkdir -p "$outdir"
    unzip "$zipfile" -d "$outdir"
    rm "$zipfile"
done
