
for dir in /media/sarlab/DATA/Bacillus_project/SRL*/; do
    if [ -d "$dir" ]; then
        shortname=$(basename "$dir")
        mkdir "$shortname"
        echo "$shortname"
    fi
done


   