project_bacillus = "/media/sarlab/DATA/Bacillus_project/"

for dir in project_bacillus/SRL*/; do
    if [ -d "$dir" ]; then
        shortname=$(basename "$dir")
        mkdir "$shortname"
        echo "$shortname"
    fi
done


   