project_bacillus = "/media/sarlab/DATA/Bacillus_project/SRL*/"

for dir in project_bacillus; do
    if [ -d "$dir" ]; then
        shortname=$(basename "$dir")
        mkdir "$shortname"
        echo "$shortname"
    fi
done


   